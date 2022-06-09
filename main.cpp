#include <iostream>
#include <fstream>
#include "solver/rocket_mpc.hpp"
#include "solver/polympc_redef.hpp" 

// #include "robot/pandaWrapper.hpp"



using namespace Eigen;
using namespace std;
using namespace std::chrono;


using admm = boxADMM<guidance_ocp::VAR_SIZE, guidance_ocp::NUM_EQ, guidance_ocp::scalar_t,
                guidance_ocp::MATRIXFMT, linear_solver_traits<guidance_ocp::MATRIXFMT>::default_solver>;

int main(int, char**) {

    // ---------- Pinocchio setup ---------- //

    // PandaWrapper myRobot;
    // Eigen::VectorXd qTarget = myRobot.inverse_kinematic(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0.5, 0., 0.5));

    
    // ---------- PolyMPC setup ---------- //

    // Creates solver
    using mpc_t = MPC<guidance_ocp, MySolver, admm>;
    mpc_t mpc;

    mpc.settings().max_iter = 20; 
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    //mpc.m_solver.settings().scaling = 10;

    // Input constraints and initialisation -------------
    const double inf = std::numeric_limits<double>::infinity();
    mpc_t::control_t lbu; 
    mpc_t::control_t ubu; 

    double max_thrust = 3000;

    lbu << 0; // lower bound on control
    ubu << max_thrust; // upper bound on control
    mpc.control_bounds(lbu, ubu);  

    // Initial control
    mpc.u_guess(ubu.replicate(mpc.ocp().NUM_NODES,1));

    // State constraints and initialisation ---------------
    double target_apogee = 2000;
    double propellant_mass = 10;
    const double eps = 1e-1;
    mpc_t::state_t lbx; 
    mpc_t::state_t ubx; 
    
    lbx << 0-eps,   0-eps,   0-eps;
    ubx <<  target_apogee+eps,      330,  propellant_mass+eps;

    mpc.state_bounds(lbx, ubx);
    
    // Final state
    mpc_t::state_t lbx_f; lbx_f << target_apogee-eps,   0-1,   0-eps; // lower bound on final state
    mpc_t::state_t ubx_f; ubx_f << target_apogee+eps,   0+1,   propellant_mass+eps; // upper bound on final state
    mpc.final_state_bounds(lbx_f, ubx_f);

    // Initial state
    mpc_t::state_t x0;
    x0 <<  0, 0, propellant_mass;
    mpc.initial_conditions(x0); 
    
    mpc.x_guess(x0.replicate(mpc.ocp().NUM_NODES,1));	
    
    // Parameters
    mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    mpc_t::parameter_t ubp; ubp << 50;   // upper bound on time
    mpc_t::parameter_t p0; p0 << 30;     // very important to set initial time estimate

    mpc.parameters_bounds(lbp, ubp);
    mpc.p_guess(p0);

    // Solve problem and save solution 
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
        std::cout << "Solution U: " << mpc.solution_u().transpose() << "\n"
        << "-------------\n";

        mpc.x_guess(mpc.solution_x());	
    }

    // Write data to txt file
    std::ofstream logFile;
    logFile.open("data/optimal_solution.txt");
    if(logFile.is_open()){

        int nPoints = 100;
        for (int iPoint = 0; iPoint<nPoints; iPoint++)
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
