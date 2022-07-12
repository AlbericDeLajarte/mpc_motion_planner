#include "motionPlanner.hpp"

MotionPlanner::MotionPlanner(): robot(), mpc(){

    Matrix<double, NDOF, 1> default_position = (robot.max_position.array() + robot.min_position.array())/2;

    set_target_state(default_position, Matrix<double, NDOF, 1>::Zero());
    set_current_state(default_position, Matrix<double, NDOF, 1>::Zero());

    // ---------- PolyMPC setup ---------- //

    // Solver settings
    mpc.settings().max_iter = 2; 
    mpc.qp_settings().max_iter = 700;
    mpc.settings().line_search_max_iter = 10;
    mpc.set_time_limits(0, 1);
    mpc.qp_settings().eps_rel = 1e-3;
    mpc.qp_settings().eps_abs = 1e-3;
    // mpc.m_solver.settings().scaling = 10;

    // Set constraints without margins (1.0) by default
    set_constraint_margins(1.0, 1.0, 1.0, 1.0, 1.0);
}

void MotionPlanner::set_target_state(Matrix<double, NDOF, 1> target_position, Matrix<double, NDOF, 1> target_velocity){

    // Update MPC constraints
    target_state.head(7) = target_position;
    target_state.tail(7) = target_velocity;

    mpc.final_state_bounds(target_state.array() - eps, target_state.array() + eps);

    // Update Ruckig constraints
    Matrix<double, 7, 1>::Map(input.target_position.data() ) = target_position;
    Matrix<double, 7, 1>::Map(input.target_velocity.data() ) = target_velocity;
    input.target_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  

}

void MotionPlanner::set_current_state(Matrix<double, NDOF, 1> current_position, Matrix<double, NDOF, 1> current_velocity){
    
    // Update MPC constraints
    current_state.head(7) = current_position;
    current_state.tail(7) = current_velocity;

    mpc.initial_conditions(current_state);


    // Update Ruckig constraints
    Matrix<double, 7, 1>::Map(input.current_position.data() ) = current_position;
    Matrix<double, 7, 1>::Map(input.current_velocity.data() ) = current_velocity;
    input.current_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  
}

void MotionPlanner::set_constraint_margins(double margin_position, double margin_velocity, double margin_acceleration, double margin_torque, double margin_jerk){

    // Save margins
    this->margin_position_ = margin_position;
    this->margin_velocity_ = margin_velocity;
    this->margin_acceleration_ = margin_acceleration;
    this->margin_torque_ = margin_torque;
    this->margin_jerk_ = margin_jerk;

    // ---------- MPC constraints ---------- //

    // State constraints ---------------
    Matrix<double, 7, 1> safety_range_position = (1-margin_position_)*(robot.max_position.array() - robot.min_position.array())/2; // needed because position bounds are not symmetric
    mpc_t::state_t lbx; lbx << robot.min_position + safety_range_position, -margin_velocity*robot.max_velocity;
    mpc_t::state_t ubx; ubx << robot.max_position - safety_range_position, margin_velocity*robot.max_velocity;
    mpc.state_bounds(lbx, ubx);

    // Input constraints -------------
    mpc.control_bounds(-margin_acceleration*robot.max_acceleration, margin_acceleration*robot.max_acceleration);  
    
    // Parameters ------------------
    mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
    mpc_t::parameter_t ubp; ubp << 10;   // upper bound on time
    mpc.parameters_bounds(lbp, ubp);

    // Non-linear torque constraints
    mpc.constraints_bounds(-margin_torque*robot.max_torque, margin_torque*robot.max_torque);


    // ---------- Ruckig constraints ---------- //
    Matrix<double, 7, 1>::Map(input.max_velocity.data() ) = margin_velocity*robot.max_velocity;
    Matrix<double, 7, 1>::Map(input.max_acceleration.data() ) = margin_acceleration*robot.max_acceleration;
    Matrix<double, 7, 1>::Map(input.max_jerk.data() ) = margin_jerk*robot.max_jerk;

}

void MotionPlanner::sample_random_state(Matrix<double, 7, 1> &random_position, Matrix<double, 7, 1> &random_velocity){
    
    Matrix<double, 7, 1> safety_range_position = (1-margin_position_)*(robot.max_position.array() - robot.min_position.array())/2;

    random_position = 0.5*(Matrix<double, 7, 1>::Random().array()*(robot.max_position - robot.min_position - 2*safety_range_position).array() 
                                                                + (robot.max_position + robot.min_position).array() );

    random_velocity = margin_velocity_*Matrix<double, 7, 1>::Random().array()*robot.max_velocity.array();
}

int MotionPlanner::check_state_in_bounds(Matrix<double, 7, 1> &position, Matrix<double, 7, 1> &velocity){

    Matrix<double, 7, 1> safety_range_position = (1-margin_position_)*(robot.max_position.array() - robot.min_position.array())/2;

    bool position_check = ( position.array() > (robot.max_position - safety_range_position).array() ).any()
                       || ( position.array() < (robot.min_position + safety_range_position).array() ).any() ;


    bool velocity_check = ( velocity.array().abs() > margin_velocity_*robot.max_velocity.array() ).any();
    
    // By default: no problem
    int check_flag = 0;

    // Position outside bounds
    if (position_check && !velocity_check) check_flag = 1;

    // Velocity outside bounds
    if (!position_check && velocity_check) check_flag = 2;

    // Position and velocity outside bounds
    if (position_check && velocity_check) check_flag = 3;

    return check_flag;
}

void MotionPlanner::warm_start_MPC(){

    // Compute Ruckig trajectory in an offline manner (outside of the control loop)
    Result result = otg.calculate(input, trajectory);

    mpc_t::traj_state_t x_guess;
    mpc_t::traj_control_t u_guess;
    mpc_t::parameter_t p0; p0 << trajectory.get_duration();
    std::array<double, NDOF> new_position, new_velocity, new_acceleration;

    auto mpc_time_grid = mpc.ocp().time_nodes;
    int i = 0;
    for(auto mpc_time : mpc_time_grid){

        trajectory.at_time(mpc_time*trajectory.get_duration(), new_position, new_velocity, new_acceleration);

        x_guess.segment(i*NDOF*2, NDOF*2) << Map<Matrix<double, NDOF, 1> >(new_position.data()),
                                                Map<Matrix<double, NDOF, 1> >(new_velocity.data());

        u_guess.segment(i*NDOF, NDOF) << Map<Matrix<double, NDOF, 1> >(new_acceleration.data());

        i++;
    } 
    // std::cout << x_guess << std::endl << x_guess.cols() << " " << x_guess.rows() << std::endl;
    // std::cout << x_guess.reshaped(14, 13)  << std::endl;
    // std::cout << u_guess.reshaped(7, 13)  << std::endl;
    mpc.x_guess(x_guess);	
    mpc.u_guess(u_guess);
    mpc.p_guess(p0); 
}

void MotionPlanner::solve_trajectory(bool use_ruckig_as_warm_start){

     // Warm start with ruckig if needed
    if (use_ruckig_as_warm_start) warm_start_MPC();

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