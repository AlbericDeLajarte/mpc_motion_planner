#include "main.hpp"


int main(int, char**) { 

    MotionPlanner planner;

    // Add margins on limits [position, velocity, acceleration, torque, jerk]
    planner.set_constraint_margins(0.7, 0.5, 0.1, 0.7, 0.01);
    

    // ---------- Compute random target state ---------- //
        
    std::srand((unsigned int) time(0));
    Matrix<double, 7, 1> target_position = planner.margin_position_*0.5*(Matrix<double, 7, 1>::Random().array()*(planner.robot.max_position - planner.robot.min_position).array() 
                                                + (planner.robot.max_position + planner.robot.min_position).array() );
    Matrix<double, 7, 1> target_velocity = planner.margin_velocity_*Matrix<double, 7, 1>::Random().array()*planner.robot.max_velocity.array();

    Matrix<double, 6, 1> task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
    

    // Check if task velocity is inside bounds, scale down otherwise
    if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) { // Linear velocity
        target_velocity *= 0.9 * planner.robot.max_linear_velocity / task_velocity.head(3).norm();
        
        std::cout << "Linear Vel: " << task_velocity.head(3).norm();
        task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
        std::cout << " corrected to : " << task_velocity.head(3).norm()  << std::endl;
    }
    if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) { // Angular velocity
        target_velocity *= 0.9 * planner.robot.max_angular_velocity / task_velocity.tail(3).norm();
        
        std::cout << "Angular Vel: " << task_velocity.tail(3).norm();
        task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
        std::cout << " corrected to : " << task_velocity.tail(3).norm() << std::endl;
    }

    std::cout << target_position.transpose() << std::endl;
    std::cout << target_velocity.transpose() << std::endl;
            
    
    // ---------- Solve trajectory ---------- //

    planner.set_target_state(target_position, target_velocity);
    planner.solve_trajectory();


    // --------- Write data to txt file --------- //
    const int nPoints = 200;

    // Storage for both trajectories
    Eigen::Matrix<double, 29, nPoints+1> ruckig_traj;
    Eigen::Matrix<double, 29, nPoints+1> polympc_traj;

    std::ofstream logFile;
    logFile.open("data/optimal_solution.txt");
    if(logFile.is_open()){

        // Log target state
        logFile << 0.0 << " " 
                << target_position.transpose() << " " 
                << target_velocity.transpose() << " "
                << Matrix<double, 1, 14>::Zero() << " " 
                << std::endl;
        
        // Get Ruckig trajectory
        Matrix<double, 7, nPoints+1> position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory;
        Matrix<double, 1, nPoints+1> time;
        planner.get_ruckig_trajectory<nPoints>(time, position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory);
        
        
        // Place data inside larger matrix for easier file writing
        ruckig_traj.row(0) = time;
        ruckig_traj.block(1, 0, 7,  nPoints+1) = position_trajectory;
        ruckig_traj.block(8, 0, 7,  nPoints+1) = velocity_trajectory;
        ruckig_traj.block(15, 0, 7,  nPoints+1) = acceleration_trajectory;
        ruckig_traj.block(22, 0, 7,  nPoints+1) = torque_trajectory;

        logFile << ruckig_traj.transpose() << std::endl;
    

        // Log Polympc trajectory
        planner.get_MPC_trajectory<nPoints>(time, position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory);
        
        // Place data inside larger matrix for easier file writing
        polympc_traj.row(0) = time;
        polympc_traj.block(1, 0, 7,  nPoints+1) = position_trajectory;
        polympc_traj.block(8, 0, 7,  nPoints+1) = velocity_trajectory;
        polympc_traj.block(15, 0, 7,  nPoints+1) = acceleration_trajectory;
        polympc_traj.block(22, 0, 7,  nPoints+1) = torque_trajectory;

        logFile << polympc_traj.transpose() << std::endl;
    }
    else std::cout << "\n !! COULD NOT OPEN FILE !!\n Data won't be saved " << std::endl;

}
