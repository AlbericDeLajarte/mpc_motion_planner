#include <iostream>
#include "rocket_mpc.hpp"
// #include "Robot_ocp.hpp"
#include "polympc_redef.hpp" 


using namespace Eigen;
using namespace std;
using namespace std::chrono;

// using box_admm_solver = boxADMM<RobotOCP::VAR_SIZE, RobotOCP::NUM_EQ, RobotOCP::scalar_t,
//                                 RobotOCP::MATRIXFMT, linear_solver_traits<RobotOCP::MATRIXFMT>::default_solver>;

// int main(void)
// {
//     /** create an MPC algorithm and set the prediction horison */
//     using mpc_t = MPC<RobotOCP, MySolver, box_admm_solver>;
//     mpc_t mpc;
//     mpc.settings().max_iter = 20;
//     mpc.settings().line_search_max_iter = 10;
//     mpc.set_time_limits(0, 2);

//     /** problem data */
//     mpc_t::static_param p; p << 2.0;          // robot wheel base
//     mpc_t::state_t x0; x0 << 0.5, 0.5, 0.5;   // initial condition
//     mpc_t::control_t lbu; lbu << -1.5, -0.75; // lower bound on control
//     mpc_t::control_t ubu; ubu <<  1.5,  0.75; // upper bound on control

//     mpc.set_static_parameters(p);
//     mpc.control_bounds(lbu, ubu);
//     mpc.initial_conditions(x0);

//     /** solve */
//     polympc::time_point start = polympc::get_time();
//     mpc.solve();
//     polympc::time_point stop = polympc::get_time();
//     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

//     /** retrieve solution and statistics */
//     std::cout << "MPC status: " << mpc.info().status.value << "\n";
//     std::cout << "Num iterations: " << mpc.info().iter << "\n";
//     std::cout << "Solve time: " << std::setprecision(9) << static_cast<double>(duration.count()) << "[mc] \n";

//     std::cout << "Solution X: " << mpc.solution_x().transpose() << "\n";
//     std::cout << "Solution U: " << mpc.solution_u().transpose() << "\n";

//     /** sample x solution at collocation points [0, 5, 10] */
//     std::cout << "x[0]: " << mpc.solution_x_at(0).transpose() << "\n";
//     std::cout << "x[5]: " << mpc.solution_x_at(5).transpose() << "\n";
//     std::cout << "x[10]: " << mpc.solution_x_at(10).transpose() << "\n";

//     std::cout << " ------------------------------------------------ \n";

//     /** sample control at collocation points */
//     std::cout << "u[0]: " << mpc.solution_u_at(0).transpose() << "\n";
//     std::cout << "u[1]: " << mpc.solution_u_at(1).transpose() << "\n";

//     std::cout << " ------------------------------------------------ \n";

//     /** sample state at time 't = [0.0, 0.5]' */
//     std::cout << "x(0.0): " << mpc.solution_x_at(0.0).transpose() << "\n";
//     std::cout << "x(0.5): " << mpc.solution_x_at(0.5).transpose() << "\n";

//     std::cout << " ------------------------------------------------ \n";

//     /**  sample control at time 't = [0.0, 0.5]' */
//     std::cout << "u(0.0): " << mpc.solution_u_at(0.0).transpose() << "\n";
//     std::cout << "u(0.5): " << mpc.solution_u_at(0.5).transpose() << "\n";

//     return EXIT_SUCCESS;
// }


using admm = boxADMM<guidance_ocp::VAR_SIZE, guidance_ocp::NUM_EQ, guidance_ocp::scalar_t,
                guidance_ocp::MATRIXFMT, linear_solver_traits<guidance_ocp::MATRIXFMT>::default_solver>;

int main(int, char**) {

    // Creates solver
    using mpc_t = MPC<guidance_ocp, MySolver, admm>;
    mpc_t mpc;

    mpc.settings().max_iter = 10; 
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    //mpc.m_solver.settings().max_iter = 1000;
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
    
    lbx << -inf, -inf, -inf,    -inf, -inf, -inf,   0-eps,   0-eps;
    ubx <<  inf,  inf, inf,      inf,  inf, inf,    propellant_mass+eps, 50;
    
    // lbx << -inf, -inf, -inf,   -inf, -inf, -inf,   -inf;
    // ubx <<  inf,  inf, inf,     inf,  inf, inf,     inf;
    mpc.state_bounds(lbx, ubx);
    
    // Final state
    mpc_t::state_t lbx_f; lbx_f << -inf, -inf, target_apogee-eps,   -inf, -inf, 0-1,   0-eps,   0-eps; // lower bound on final state
    mpc_t::state_t ubx_f; ubx_f <<  inf,  inf, target_apogee+eps,    inf,  inf, 0+1,   propellant_mass+eps, 50; // upper bound on final state
    mpc.final_state_bounds(lbx_f, ubx_f);

    // Initial state
    mpc_t::state_t x0_inf, x0_sup;
    x0_inf <<   0, 0, 0,
                0, 0, 0,
                propellant_mass,
                0;
    x0_sup = x0_inf;
    x0_sup(7) = 50;
    mpc.initial_conditions(x0_inf, x0_sup); 
    
    mpc.x_guess(x0_sup.replicate(mpc.ocp().NUM_NODES,1));	
    
    // Parameters
    // mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    // mpc_t::parameter_t ubp; ubp << 50;   // upper bound on time
    // mpc_t::parameter_t p0; p0 << 30;     // very important to set initial time estimate

    // mpc.parameters_bounds(lbp, ubp);
    // mpc.p_guess(p0);

    mpc.ocp().cos_vertical_angle = 1;

    // Solve problem and save solution 
    for(int i=0; i<10; i++){
        auto start = system_clock::now();

        mpc.solve(); 
        
        auto stop = system_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        double dT = duration.count()*1e-6;
       
        /** retrieve solution and statistics */
        std::cout << "MPC status: " << mpc.info().status.value << "\n";
        std::cout << "Num iterations: " << mpc.info().iter << "\n";
        std::cout << "Solve time: " << dT << " [s] \n";

        // std::cout << "Final time: " << mpc.solution_p().transpose() << endl;

        std::cout << "Solution X: \n" << mpc.solution_x().reshaped(8, 11).transpose() << "\n";
        std::cout << "Solution U: " << mpc.solution_u().transpose() << "\n"
        << "-------------\n";

        // mpc.x_guess(mpc.solution_x());	
    }

}
