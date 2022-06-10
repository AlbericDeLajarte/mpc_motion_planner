#include <iostream>
#include <fstream>
#include "mpc_solver/robot_ocp.hpp"
#include "mpc_solver/polympc_redef.hpp" 

#include "robot/pandaWrapper.hpp"



using namespace Eigen;
using namespace std;
using namespace std::chrono;


using admm = boxADMM<minTime_ocp::VAR_SIZE, minTime_ocp::NUM_EQ, minTime_ocp::scalar_t,
                minTime_ocp::MATRIXFMT, linear_solver_traits<minTime_ocp::MATRIXFMT>::default_solver>;

int main(int, char**) {

    // ---------- PolyMPC setup ---------- //

    // Creates solver
    using mpc_t = MPC<minTime_ocp, MySolver, admm>;
    mpc_t mpc;

    mpc.settings().max_iter = 10; 
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    // mpc.m_solver.settings().scaling = 10;

    // State constraints and initialisation ---------------
    mpc_t::state_t lbx; 
    mpc_t::state_t ubx; 
    
    // Limits from https://frankaemika.github.io/docs/control_parameters.html
    lbx << -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973,
            -2.1750, -2.1750, -2.1750, -2.1750, -2.6100, -2.6100, -2.6100;
    ubx <<  2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973,
            2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100;
    mpc.state_bounds(lbx, ubx);

    // Input constraints and initialisation -------------
    const double inf = std::numeric_limits<double>::infinity();
    mpc_t::control_t max_acceleration; 

    max_acceleration << 15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0; // acceleration limit
    mpc.control_bounds(-max_acceleration, max_acceleration);  
    // mpc.u_guess(ubu.replicate(mpc.ocp().NUM_NODES,1));

    
    // Parameters ------------------
    mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    mpc_t::parameter_t ubp; ubp << 10;   // upper bound on time
    mpc_t::parameter_t p0; p0 << 3;     // very important to set initial time estimate

    mpc.parameters_bounds(lbp, ubp);
    mpc.p_guess(p0);    



    // ---------- Pinocchio setup ---------- //

    PandaWrapper myRobot;

    bool feasibleTarget = false;
    Matrix<double, 7, 1> qTarget;
    mpc_t::state_t final_state; 
    
    // Search over target configuration until one is inside joint limits
    while (feasibleTarget == false){
        qTarget = myRobot.inverse_kinematic(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0.5, 0., 0.5));
        cout << qTarget.transpose() << endl;

        // Final state
        final_state.head(7) = qTarget;
        final_state.tail(7) = Eigen::Matrix<double, 7, 1>::Zero();

        // Check state constraints violation
        if( (final_state.array() < ubx.array()).all() && (final_state.array() > lbx.array()).all() ){
            cout << "OK" << endl;
            feasibleTarget = true;
        }
         
        else cout << "NOT OK" << endl;
    }

    
    // Constraint initial and final states ---------------
    const double eps = 1e-2;
    mpc.final_state_bounds(final_state.array() - eps, final_state.array() + eps);

    mpc_t::state_t x0; x0 << 0.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    mpc.initial_conditions(x0);

    mpc.x_guess(x0.replicate(mpc.ocp().NUM_NODES,1));	
    

    // Solve problem and print solution 
    for(int i=0; i<5; i++){
        auto start = system_clock::now();

        mpc.solve(); 
        
        auto stop = system_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        double dT = duration.count()*1e-3;
       
        /** retrieve solution and statistics */
        std::cout << "MPC status: " << mpc.info().status.value << "\n";
        std::cout << "Num iterations: " << mpc.info().iter << "\n";
        std::cout << "Solve time: " << dT << " [ms] \n";

        std::cout << "Final time: " << mpc.solution_p().transpose() << endl;

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
        for (int iPoint = 0; iPoint<nPoints; iPoint++)
        {
            double time = 1.0/nPoints * iPoint;
            
            // Log optimal trajectory
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
