#include "main.hpp"
#include "pinocchio/algorithm/rnea.hpp"

using namespace Eigen;

using namespace ruckig;

using admm = boxADMM<minTime_ocp::VAR_SIZE, minTime_ocp::NUM_EQ + minTime_ocp::NUM_INEQ, minTime_ocp::scalar_t,
                minTime_ocp::MATRIXFMT, linear_solver_traits<minTime_ocp::MATRIXFMT>::default_solver>;


int main(int, char**) { 

    PandaWrapper robot;
    MotionPlanner planner;

    // Reducing a lot acceleration limits to force it to use longer path 
    // and generate potential joint limit violations
    for(int i = 0; i< 7; i++){
        robot.max_acceleration.at(i) /= 1;
    }

    // ---------- PolyMPC setup ---------- //

    // Creates solver
    using mpc_t = MPC<minTime_ocp, MySolver, admm>;
    mpc_t mpc;

    // mpc.ocp().model = robot.model;
    // mpc.ocp().data = robot.data;

    mpc.settings().max_iter = 3; 
    mpc.qp_settings().max_iter = 100;
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    mpc.qp_settings().eps_rel = 1e-3;
    mpc.qp_settings().eps_abs = 1e-3;
    // mpc.m_solver.settings().scaling = 10;

    std::cout << "rel: " << mpc.qp_settings().eps_rel << " abs: " <<  mpc.qp_settings().eps_abs << std::endl;

    // State constraints and initialisation ---------------
    mpc_t::state_t lbx; 
    mpc_t::state_t ubx; 
    
    // Limits from https://frankaemika.github.io/docs/control_parameters.html
    lbx << Map<Matrix<double, 7, 1> >(robot.min_position.data()),
          -Map<Matrix<double, 7, 1> >(robot.max_velocity.data());
    ubx << Map<Matrix<double, 7, 1> >(robot.max_position.data()),
           Map<Matrix<double, 7, 1> >(robot.max_velocity.data());
    mpc.state_bounds(lbx, ubx);

    // Input constraints and initialisation -------------
    const double inf = std::numeric_limits<double>::infinity();
    mpc_t::control_t max_input; 

    max_input = Map<Matrix<double, 7, 1> >(robot.max_acceleration.data()); // acceleration limit
    mpc.control_bounds(-max_input, max_input);  
    
    // Parameters ------------------
    mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    mpc_t::parameter_t ubp; ubp << 10;   // upper bound on time

    mpc.parameters_bounds(lbp, ubp);

    // Non-linear constraints
    mpc_t::constraint_t ubg, lbg;
    lbg << -inf;
    ubg << inf;

    mpc.constraints_bounds(lbg, ubg);
    



    // ---------- Pinocchio setup ---------- //

    bool feasibleTarget = false;
    Matrix<double, 7, 1> qTarget;
    mpc_t::state_t final_state; 
    
    // Search over target configuration until one is inside joint limits
    int nTry = 0;
    while (feasibleTarget == false && nTry < 100){
        qTarget = robot.inverse_kinematic(Matrix3d::Identity(), Vector3d(0.4, 0., 0.5));

        // Final state
        final_state.head(7) = qTarget;

        // Check state constraints violation
        if( (final_state.array() < ubx.array()).all() && (final_state.array() > lbx.array()).all() ){
            std::cout << "Solved in " << nTry << " trials. Initial state: " << std::endl;
            feasibleTarget = true;
        }
        nTry ++;
    }

    // Compute desired final joint speed from cartesian [linear, angular] speed
    final_state.tail(7) = robot.inverse_velocities(qTarget, Vector3d(0.5, 0., 0.3), Vector3d(0.0, 0.0, 0.0));

    std::cout << final_state.reshaped(7, 2).transpose() << std::endl;

    std::cout <<  pinocchio::rnea(robot.model, robot.data, qTarget, qTarget, qTarget);
    

    // ---------- Ruckig setup ---------- //
    // Create input parameters
    InputParameter<NDOF> input;
    input.current_position = planner.init_position;
    input.current_velocity = planner.init_velocity;
    input.current_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Matrix<double, 7, 1>::Map(input.target_position.data() ) = final_state.head(7);
    Matrix<double, 7, 1>::Map(input.target_velocity.data() ) = final_state.tail(7);
    input.target_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    input.max_velocity = robot.max_velocity;
    input.max_acceleration = robot.max_acceleration;
    input.max_jerk = robot.max_jerk;

    // We don't need to pass the control rate (cycle time) when using only offline features
    Ruckig<NDOF> otg;
    Trajectory<NDOF> trajectory;

    // Calculate the trajectory in an offline manner (outside of the control loop)
    Result result = otg.calculate(input, trajectory);

    // Get duration of the trajectory
    std::cout << "Ruckig trajectory duration: " << trajectory.get_duration() << " [s]. \n\n";

    
    // ---------- SOLVE POLYMPC ---------- //

    // Constraint initial and final state ---------------
    const double eps = 1e-2;
    mpc.final_state_bounds(final_state.array() - eps, final_state.array() + eps);

    mpc_t::state_t x0; x0 << Map<Matrix<double, 7, 1> >(planner.init_position.data()),
                             Map<Matrix<double, 7, 1> >(planner.init_velocity.data());
    mpc.initial_conditions(x0);

    // Warm start solver
    mpc_t::traj_state_t x_guess;
    mpc_t::traj_control_t u_guess;
    mpc_t::parameter_t p0; p0 << trajectory.get_duration();
    std::array<double, NDOF> new_position, new_velocity, new_acceleration;

    auto mpc_time_grid = mpc.ocp().time_nodes;
    int i = 0;
    for(auto mpc_time : mpc_time_grid){

        trajectory.at_time(mpc_time*trajectory.get_duration(), new_position, new_velocity, new_acceleration);

        // std::cout << mpc_time*trajectory.get_duration() << std::endl;

        x_guess.segment(i*NDOF*2, NDOF*2) << Map<Matrix<double, 7, 1> >(new_position.data()),
                                             Map<Matrix<double, 7, 1> >(new_velocity.data());

        u_guess.segment(i*NDOF, NDOF) << Map<Matrix<double, 7, 1> >(new_acceleration.data());

        i++;
    } 
    // std::cout << x_guess << std::endl << x_guess.cols() << " " << x_guess.rows() << std::endl;
    // std::cout << x_guess.reshaped(14, 13)  << std::endl;
    // std::cout << u_guess.reshaped(7, 13)  << std::endl;
    mpc.x_guess(x_guess);	
    mpc.u_guess(u_guess);
    mpc.p_guess(p0); 
     

    // Solve problem and print solution 
    std::cout << " ---------- SOLVING MPC ----------" << std::endl;
    for(int i=0; i<2; i++){
        auto start = std::chrono::system_clock::now();

        mpc.solve(); 
        
        auto stop = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        double dT = duration.count()*1e-3;
       
        /** retrieve solution and statistics */
        std::cout << "MPC status: " << mpc.info().status.value << "\n";
        std::cout << "Num iterations: " << mpc.info().iter << "\n";
        std::cout << "Solve time: " << dT << " [ms] \n";

        std::cout << "Final time: " << mpc.solution_p().transpose() << std::endl;

        // std::cout << "Solution X: \n" << mpc.solution_x().reshaped(3, 6).transpose() << "\n";
        // std::cout << "Solution U: " << mpc.solution_u().transpose() << "\n"
        std::cout << "-------------\n";

        mpc.x_guess(mpc.solution_x());	
        mpc.u_guess(mpc.solution_u());
        mpc.p_guess(mpc.solution_p());
    }

    // Write data to txt file
    std::ofstream logFile;
    logFile.open("data/optimal_solution.txt");
    if(logFile.is_open()){

        // Log target state
        logFile << 0.0 << " " 
                << final_state.transpose() << " " 
                << Matrix<double, 1, 7>::Zero() << " " 
                << std::endl;

        int nPoints = 100;

        // Log Ruckig trajectory
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            double time = trajectory.get_duration()/nPoints * iPoint;
            
            trajectory.at_time(time, new_position, new_velocity, new_acceleration);
            
            logFile << time << " " 
                    << Map<Matrix<double, 1, 7> >(new_position.data()) << " " 
                    << Map<Matrix<double, 1, 7> >(new_velocity.data()) << " " 
                    << Map<Matrix<double, 1, 7> >(new_acceleration.data()) << " " 
                    << std::endl;
        }

        // Log Polympc trajectory
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            double time = 1.0/nPoints * iPoint;
            
            logFile << time*mpc.solution_p()[0] << " " 
                    << mpc.solution_x_at(time).transpose() << " " 
                    << mpc.solution_u_at(time).transpose() << " " 
                    << std::endl;
        }
    }
    else {
        std::cout << "\n !! COULD NOT OPEN FILE !!\n Data won't be saved " << std::endl;
    }

}
