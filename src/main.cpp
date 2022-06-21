#include "main.hpp"
#include "pinocchio/algorithm/rnea.hpp"

using namespace Eigen;

using namespace ruckig;

using admm = boxADMM<minTime_ocp::VAR_SIZE, minTime_ocp::NUM_EQ + minTime_ocp::NUM_INEQ, minTime_ocp::scalar_t,
                minTime_ocp::MATRIXFMT, linear_solver_traits<minTime_ocp::MATRIXFMT>::default_solver>;


int main(int, char**) { 

    PandaWrapper robot;
    MotionPlanner planner;

    robot.min_position = robot.min_position.array()+0.3;
    robot.max_position = robot.max_position.array()-0.3;

    robot.max_velocity *= 0.3;
    robot.max_acceleration *= 0.7;
    robot.max_jerk *= 0.01;
    robot.max_torque *= 1;

    robot.max_torqueDot *= 1;
    // ---------- PolyMPC setup ---------- //

    // Creates solver
    using mpc_t = MPC<minTime_ocp, MySolver, admm>;
    mpc_t mpc;

    mpc.settings().max_iter = 2; 
    mpc.qp_settings().max_iter = 700;
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    mpc.qp_settings().eps_rel = 1e-3;
    mpc.qp_settings().eps_abs = 1e-3;
    // mpc.m_solver.settings().scaling = 10;

    // State constraints and initialisation ---------------
    mpc_t::state_t lbx; 
    mpc_t::state_t ubx; 
    
    // Limits from https://frankaemika.github.io/docs/control_parameters.html
    lbx << robot.min_position, -robot.max_velocity;
    ubx << robot.max_position, robot.max_velocity;
    mpc.state_bounds(lbx, ubx);

    // Input constraints and initialisation -------------
    const double inf = std::numeric_limits<double>::infinity();
    mpc_t::control_t max_input; 

    max_input = robot.max_acceleration; // acceleration limit
    mpc.control_bounds(-max_input, max_input);  
    
    // Parameters ------------------
    mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    mpc_t::parameter_t ubp; ubp << 10;   // upper bound on time

    mpc.parameters_bounds(lbp, ubp);

    // Non-linear torque constraints
    mpc_t::constraint_t ubg, lbg;
    lbg = -robot.max_torque;
    ubg = robot.max_torque;
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
        // final_state.head(7) = qTarget;

        // Check state constraints violation
        if( (final_state.array() < ubx.array()).all() && (final_state.array() > lbx.array()).all() ){
            std::cout << "Solved in " << nTry << " trials. Initial state: " << std::endl;
            feasibleTarget = true;
        }
        nTry ++;
    }

    // qTarget << 3.14*0.5, -0.936495238095238, -2.3686368421052633, -0.9453, 2.368636842105263, 2.9079545454545457, 0.4281842865538632;
    // final_state.head(7) = qTarget;
    // final_state.tail(7) << 1.80852545,  1.78654623,  1.53257108,  1.82629974, -0.03752599, -1.46085578 , 0.0;
    // final_state.tail(7) *= 0.8;

    for(int iter =0; iter<1000; iter++){

    final_state.head(7) = 0.5*(Matrix<double, 7, 1>::Random().array()*(robot.max_position-robot.min_position).array() + (robot.max_position+robot.min_position).array() );
    final_state.tail(7) = Matrix<double, 7, 1>::Random().array()*robot.max_velocity.array();
    // // Compute desired final joint speed from cartesian [linear, angular] speed
    // final_state.tail(7) = robot.inverse_velocities(qTarget, Vector3d(0.5, 0., 0.3), Vector3d(0.0, 0.0, 0.0));

    std::cout << final_state.reshaped(7, 2).transpose() << std::endl;    


    // ---------- Ruckig setup ---------- //
    // Create input parameters
    InputParameter<NDOF> input;

    Matrix<double, 7, 1>::Map(input.current_position.data() ) = planner.init_position;
    Matrix<double, 7, 1>::Map(input.current_velocity.data() ) = planner.init_velocity;
    input.current_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Matrix<double, 7, 1>::Map(input.target_position.data() ) = final_state.head(7);
    Matrix<double, 7, 1>::Map(input.target_velocity.data() ) = final_state.tail(7);
    input.target_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Matrix<double, 7, 1>::Map(input.max_velocity.data() ) = robot.max_velocity;
    Matrix<double, 7, 1>::Map(input.max_acceleration.data() ) = robot.max_acceleration;
    Matrix<double, 7, 1>::Map(input.max_jerk.data() ) = robot.max_jerk;

    // We don't need to pass the control rate (cycle time) when using only offline features
    Ruckig<NDOF> otg;
    Trajectory<NDOF> trajectory;

    // Calculate the trajectory in an offline manner (outside of the control loop)
    Result result = otg.calculate(input, trajectory);

    // Get duration of the trajectory
    // std::cout << "Ruckig trajectory duration: " << trajectory.get_duration() << " [s]. \n\n";

    
    // ---------- SOLVE POLYMPC ---------- //

    // Constraint initial and final state ---------------
    const double eps = 1e-2;
    mpc.final_state_bounds(final_state.array() - eps, final_state.array() + eps);

    mpc_t::state_t x0; x0 << planner.init_position, planner.init_velocity;
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
    // std::cout << " ---------- SOLVING MPC ----------" << std::endl;
    for(int i=0; i<1; i++){
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

    // --------- Write data to txt file --------- //
    const int nPoints = 200;
    Eigen::Matrix<double, 28, nPoints+1> ruckig_traj;
    Eigen::Matrix<double, 28, nPoints+1> polympc_traj;

    int linear_vel_flag_rk {1}, angular_vel_flag_rk {1}, jerk_flag_rk {1}, torqueDot_flag_rk {1};
    int linear_vel_flag_mpc {1}, angular_vel_flag_mpc {1}, jerk_flag_mpc {1}, torqueDot_flag_mpc {1};

    // Save trajectories
    std::ofstream logFile;
    logFile.open("data/optimal_solution.txt");
    if(logFile.is_open()){

        // Log target state
        logFile << 0.0 << " " 
                << final_state.transpose() << " " 
                << Matrix<double, 1, 14>::Zero() << " " 
                << std::endl;

        // Log Ruckig trajectory
        double dT = trajectory.get_duration()/nPoints;
        
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            double time = dT * iPoint;
            
            trajectory.at_time(time, new_position, new_velocity, new_acceleration);

            ruckig_traj.col(iPoint).head(21) << Map<Matrix<double, 7, 1> >(new_position.data()),
                                                Map<Matrix<double, 7, 1> >(new_velocity.data()),
                                                Map<Matrix<double, 7, 1> >(new_acceleration.data());

            ruckig_traj.col(iPoint).tail(7) = 
            pinocchio::rnea(robot.model, robot.data, ruckig_traj.col(iPoint).head(7),
                                                     ruckig_traj.col(iPoint).segment(7, 7),
                                                     ruckig_traj.col(iPoint).segment(14, 7));

            // Check jerk and torque variation
            if (iPoint >= 1){
                Matrix<double, 7, 1> jerk = (ruckig_traj.col(iPoint).segment(14, 7)-ruckig_traj.col(iPoint-1).segment(14, 7)) / dT;
                Matrix<double, 7, 1> torqueDot = (ruckig_traj.col(iPoint).tail(7)-ruckig_traj.col(iPoint-1).tail(7)) / dT;

                if( (jerk.array().abs() > 10*robot.max_jerk.array()).any() ){ // Using 10% margin because RK is driving at the limit
                    jerk_flag_rk = 0;
                    // std::cout << "RK: Jerk limit of: " << jerk.transpose() << " at time: " << time << std::endl;
                } 
                if( (torqueDot.array().abs() > robot.max_torqueDot).any() ){
                    torqueDot_flag_rk = 0;
                    // std::cout << "RK: TorqueDot limit of: " << torqueDot.transpose() << " at time: " << time << std::endl;
                } 
            }

            // Check cartesian velocity
            pinocchio::Data::Matrix6x J(6,7); J.setZero();
            pinocchio::computeJointJacobian(robot.model, robot.data, ruckig_traj.col(iPoint).head(7), 7, J);
            Matrix<double, 6, 1> task_velocity = J*ruckig_traj.col(iPoint).segment(7, 7);

            if(task_velocity.head(3).norm() > robot.max_linear_velocity) {
                linear_vel_flag_rk = 0;
                // std::cout << "RK: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << time << std::endl;
            }
            if(task_velocity.tail(3).norm() > robot.max_angular_velocity) {
                angular_vel_flag_rk = 0;
                // std::cout << "RK: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << time << std::endl;
            }
            
            logFile << time << " " 
                    << ruckig_traj.col(iPoint).transpose() 
                    << std::endl;
        }

        // Log Polympc trajectory
        dT = mpc.solution_p()[0]/nPoints;
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            double time = 1.0/nPoints * iPoint;

            polympc_traj.col(iPoint).head(21) << mpc.solution_x_at(time), mpc.solution_u_at(time);

            polympc_traj.col(iPoint).tail(7) = 
            pinocchio::rnea(robot.model, robot.data, polympc_traj.col(iPoint).head(7),
                                                     polympc_traj.col(iPoint).segment(7, 7),
                                                     polympc_traj.col(iPoint).segment(14, 7));
            
            // Check jerk and torque variation
            if (iPoint >= 1){
                Matrix<double, 7, 1> jerk = (polympc_traj.col(iPoint).segment(14, 7)-polympc_traj.col(iPoint-1).segment(14, 7)) / dT;
                Matrix<double, 7, 1> torqueDot = (polympc_traj.col(iPoint).tail(7)-polympc_traj.col(iPoint-1).tail(7)) / dT;

                if( (jerk.array().abs() > 10*robot.max_jerk.array()).any() ){
                    jerk_flag_mpc = 0;
                    // std::cout << "MPC: Jerk limit of: " << jerk.transpose() << " at time: " << time << std::endl;
                } 
                if( (torqueDot.array().abs() > robot.max_torqueDot).any() ){
                    torqueDot_flag_mpc = 0;
                    std::cout << "MPC: TorqueDot limit of: " << torqueDot.transpose() << " at time: " << time << std::endl;
                } 
            }

            // Check cartesian velocity
            pinocchio::Data::Matrix6x J(6,7); J.setZero();
            pinocchio::computeJointJacobian(robot.model, robot.data, polympc_traj.col(iPoint).head(7), 7, J);
            Matrix<double, 6, 1> task_velocity = J*polympc_traj.col(iPoint).segment(7, 7);

            if(task_velocity.head(3).norm() > robot.max_linear_velocity) {
                linear_vel_flag_mpc = 0;
                // std::cout << "MPC: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << time << std::endl;
            }
            if(task_velocity.tail(3).norm() > robot.max_angular_velocity) {
                angular_vel_flag_mpc = 0;
                // std::cout << "MPC: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << time << std::endl;
            }

            logFile << time*mpc.solution_p()[0] << " " 
                    << polympc_traj.col(iPoint).transpose() 
                    << std::endl;
        }
    }
    else {
        std::cout << "\n !! COULD NOT OPEN FILE !!\n Data won't be saved " << std::endl;
    }

    // Save trajectory performance for benchmark

    std::ofstream benchFile;
    benchFile.open("data/benchmark_data.txt", std::ios_base::app);
    if(benchFile.is_open()){

        // Log Ruckig data  
        trajectory.at_time(trajectory.get_duration(), new_position, new_velocity, new_acceleration);

        // Write extremum of both trajectories
        for(int i=0; i< 28; i++) benchFile << ruckig_traj.row(i).minCoeff() << " ";
        for(int i=0; i< 28; i++) benchFile << ruckig_traj.row(i).maxCoeff() << " ";
        for(int i=0; i< 28; i++) benchFile << polympc_traj.row(i).minCoeff() << " ";
        for(int i=0; i< 28; i++) benchFile << polympc_traj.row(i).maxCoeff() << " ";

        // Write final state
        benchFile << (ruckig_traj.col(nPoints).head(14)-final_state).transpose() << " "
                  << (polympc_traj.col(nPoints).head(14)-final_state).transpose() << " ";
            
        benchFile << jerk_flag_rk << " " 
                  << torqueDot_flag_rk << " " 
                  << linear_vel_flag_rk << " " 
                  << angular_vel_flag_rk << " " 

                  << jerk_flag_mpc << " " 
                  << torqueDot_flag_mpc << " "
                  << linear_vel_flag_mpc << " " 
                  << angular_vel_flag_mpc << std::endl;
    }

    }

}
