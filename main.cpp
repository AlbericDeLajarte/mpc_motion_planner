#include <iostream>
#include <fstream>
#include "mpc_solver/robot_ocp.hpp"
#include "mpc_solver/polympc_redef.hpp" 

#include "robot/pandaWrapper.hpp"



using namespace Eigen;
using namespace std;
using namespace std::chrono;


using admm = boxADMM<guidance_ocp::VAR_SIZE, guidance_ocp::NUM_EQ, guidance_ocp::scalar_t,
                guidance_ocp::MATRIXFMT, linear_solver_traits<guidance_ocp::MATRIXFMT>::default_solver>;

int main(int, char**) {

    // ---------- Pinocchio setup ---------- //

    // PandaWrapper myRobot;
    // Eigen::VectorXd qTarget = myRobot.inverse_kinematic(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0.5, 0., 0.5));
    Eigen::Matrix<float, 7, 1> qTarget; qTarget << 0.9066, 1.16943, 1.44188, -1.31434, 1.20663, 1.786, 0.473774;
    
    // ---------- PolyMPC setup ---------- //

    // Creates solver
    using mpc_t = MPC<guidance_ocp, MySolver, admm>;
    mpc_t mpc;

    double final_time = 3;

    mpc.ocp().Q.setIdentity();
    mpc.ocp().R.setZero();
    mpc.ocp().QN.setIdentity();

    // mpc.ocp().x_ref

    mpc.settings().max_iter = 10; 
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, final_time);
    // mpc.m_solver.settings().scaling = 10;

    // Input constraints and initialisation -------------
    const float inf = std::numeric_limits<float>::infinity();
    mpc_t::control_t max_acceleration; 

    max_acceleration << 15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0; // acceleration limit
    mpc.control_bounds(-max_acceleration, max_acceleration);  

    // Initial control
    // mpc.u_guess(ubu.replicate(mpc.ocp().NUM_NODES,1));

    // State constraints and initialisation ---------------
    const float eps = 1e-2;
    mpc_t::state_t lbx; 
    mpc_t::state_t ubx; 
    
    // Position and velocity limits from https://frankaemika.github.io/docs/control_parameters.html
    lbx << -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973,
            -2.1750, -2.1750, -2.1750, -2.1750, -2.6100, -2.6100, -2.6100;
    ubx <<  2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973,
            2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100;
 
    mpc.state_bounds(lbx, ubx);
    
    // Final state
    // mpc_t::state_t final_state; 
    // final_state.head(7) = qTarget;
    // final_state.tail(7) = Eigen::Matrix<float, 7, 1>::Zero();

    mpc.ocp().x_ref.head(7) = qTarget;
    mpc.ocp().x_ref.tail(7) = Eigen::Matrix<float, 7, 1>::Zero();

    // mpc.final_state_bounds(final_state.array() - eps, final_state.array() + eps);

    // Initial state
    mpc_t::state_t x0; x0 << 0.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    mpc.initial_conditions(x0);

    mpc.x_guess(x0.replicate(mpc.ocp().NUM_NODES,1));	
    
    // Parameters
    // mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    // mpc_t::parameter_t ubp; ubp << 20;   // upper bound on time
    // mpc_t::parameter_t p0; p0 << 5;     // very important to set initial time estimate

    // mpc.parameters_bounds(lbp, ubp);
    // mpc.p_guess(p0);

    // Solve problem and save solution 
    for(int i=0; i<5; i++){
        auto start = system_clock::now();

        mpc.solve(); 
        
        auto stop = system_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        float dT = duration.count()*1e-3;
       
        /** retrieve solution and statistics */
        std::cout << "MPC status: " << mpc.info().status.value << "\n";
        std::cout << "Num iterations: " << mpc.info().iter << "\n";
        std::cout << "Solve time: " << dT << " [ms] \n";

        std::cout << "Final time: " << mpc.solution_p().transpose() << endl;

        // std::cout << "Solution X: \n" << mpc.solution_x().reshaped(3, 6).transpose() << "\n";
        std::cout << "Solution U: " << mpc.solution_u().transpose() << "\n"
        << "-------------\n";

        mpc.x_guess(mpc.solution_x());	
        mpc.u_guess(mpc.solution_u());
    }

    // Write data to txt file
    std::ofstream logFile;
    logFile.open("data/optimal_solution.txt");
    if(logFile.is_open()){

        int nPoints = 100;
        for (int iPoint = 0; iPoint<nPoints; iPoint++)
        {
            float time = final_time/nPoints * iPoint;

            logFile << time << " " 
                    << mpc.solution_x_at(time).transpose() << " " 
                    << mpc.solution_u_at(time).transpose() << " " 
                    << std::endl;
        }
    }
    else {
        std::cout << "\n !! COULD NOT OPEN FILE !!\n Data won't be saved " << std::endl;
    }

}
