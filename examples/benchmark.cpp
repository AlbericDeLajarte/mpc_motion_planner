#include "motionPlanner.hpp"


int main(int, char**) { 

    MotionPlanner planner;

    // Add margins on limits [position, velocity, acceleration, torque, jerk]
    planner.set_constraint_margins(0.9, 0.9, 0.5, 0.9, 0.1);
    

    // ---------- Compute random target state ---------- //
    std::srand((unsigned int) time(0));
    Matrix<double, 7, 1> target_position, target_velocity;

    for(int iter =0; iter<1000; iter++){   

        // Sample random state. Redefine target velocity to have zero angular speed
        planner.sample_random_state(target_position, target_velocity);
        target_velocity = planner.robot.inverse_velocities(target_position, Matrix<double, 3, 1>::Random()*planner.robot.max_linear_velocity, Vector3d(0.0, 0.0, 0.0));

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
        // Check if joint velocity is inside bounds, scale down otherwise
        if( (target_velocity.array().abs() > planner.margin_velocity_*planner.robot.max_velocity.array()).any() ){ 
            target_velocity /= 1.1*(target_velocity.array().abs()/(planner.margin_velocity_*planner.robot.max_velocity.array())).maxCoeff();
        }
     
        std::cout << target_position.transpose() << std::endl;
        std::cout << target_velocity.transpose() << std::endl;
            
        
        // ---------- Solve trajectory ---------- //

        planner.set_target_state(target_position, target_velocity);
        
        bool use_ruckig_as_warm_start = true;
        planner.solve_trajectory(use_ruckig_as_warm_start);
        

        // --------- Check trajectories and save results --------- //

        const int nPoints = 200;

        // Storage for both trajectories
        Eigen::Matrix<double, 29, nPoints+1> ruckig_traj;
        Eigen::Matrix<double, 29, nPoints+1> polympc_traj;

        Matrix<double, 7, nPoints+1> position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory;
        Matrix<double, 1, nPoints+1> time;

        int linear_vel_flag_rk {1}, angular_vel_flag_rk {1}, jerk_flag_rk {1}, torqueDot_flag_rk {1}, collision_flag_rk {1};
        int linear_vel_flag_mpc {1}, angular_vel_flag_mpc {1}, jerk_flag_mpc {1}, torqueDot_flag_mpc {1}, collision_flag_mpc {1};
        
        // --------- RUCKIG --------- //
        // Get trajectory
        planner.get_ruckig_trajectory<nPoints>(time, position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory);

        // Place data inside larger matrix for easier check
        ruckig_traj.row(0) = time;
        ruckig_traj.block(1, 0, 7,  nPoints+1) = position_trajectory;
        ruckig_traj.block(8, 0, 7,  nPoints+1) = velocity_trajectory;
        ruckig_traj.block(15, 0, 7,  nPoints+1) = acceleration_trajectory;
        ruckig_traj.block(22, 0, 7,  nPoints+1) = torque_trajectory;
        
        // Compute jerk and cartesian velocity
        double dT = time(nPoints)/nPoints;
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            // Check jerk 
            if (iPoint >= 1){
                Matrix<double, 7, 1> jerk = (acceleration_trajectory.col(iPoint)-acceleration_trajectory.col(iPoint-1)) / dT;

                if( (jerk.array().abs() > 10*planner.robot.max_jerk.array()).any() ){ // Using 10% margin because RK is driving at the limit
                    jerk_flag_rk = 0;
                    // std::cout << "RK: Jerk limit of: " << jerk.transpose() << " at time: " << dT*iPoint << std::endl;
                } 
            }

            // Check cartesian velocity
            task_velocity = planner.robot.forward_velocities(position_trajectory.col(iPoint), velocity_trajectory.col(iPoint));

            if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) {
                linear_vel_flag_rk = 0;
                // std::cout << "RK: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << dT*iPoint << std::endl;
            }
            if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) {
                angular_vel_flag_rk = 0;
                // std::cout << "RK: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << dT*iPoint << std::endl;
            }

            // Check table collision
            pinocchio::forwardKinematics(planner.robot.model, planner.robot.data, position_trajectory.col(iPoint));
            if(planner.robot.data.oMi[7].translation()[2]< 0.16) {
                collision_flag_rk = 0;
                std::cout << "RK Table collision: " << planner.robot.data.oMi[7].translation()[2] << std::endl;
            }
        }

        // --------- MPC --------- //
        // Check trajectory
        planner.get_MPC_trajectory<nPoints>(time, position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory);

        // Place data inside larger matrix for easier check
        polympc_traj.row(0) = time;
        polympc_traj.block(1, 0, 7,  nPoints+1) = position_trajectory;
        polympc_traj.block(8, 0, 7,  nPoints+1) = velocity_trajectory;
        polympc_traj.block(15, 0, 7,  nPoints+1) = acceleration_trajectory;
        polympc_traj.block(22, 0, 7,  nPoints+1) = torque_trajectory;
        
        // Compute jerk and cartesian velocity
        dT = time(nPoints)/nPoints;
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {            
            // Check jerk 
            if (iPoint >= 1){
                Matrix<double, 7, 1> jerk = (acceleration_trajectory.col(iPoint)-acceleration_trajectory.col(iPoint-1)) / dT;

                if( (jerk.array().abs() > 10*planner.robot.max_jerk.array()).any() ){
                    jerk_flag_mpc = 0;
                    // std::cout << "MPC: Jerk limit of: " << jerk.transpose() << " at time: " << dT*iPoint << std::endl;
                } 
            }

            // Check cartesian velocity
            task_velocity = planner.robot.forward_velocities(position_trajectory.col(iPoint), velocity_trajectory.col(iPoint));

            if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) {
                linear_vel_flag_mpc = 0;
                // std::cout << "MPC: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << dT*iPoint << std::endl;
            }
            if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) {
                angular_vel_flag_mpc = 0;
                // std::cout << "MPC: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << dT*iPoint << std::endl;
            }

            // Check table collision
            pinocchio::forwardKinematics(planner.robot.model, planner.robot.data, position_trajectory.col(iPoint));
            if(planner.robot.data.oMi[7].translation()[2]< 0.16) {
                collision_flag_mpc = 0;
                std::cout << "MPC Table collision: " << planner.robot.data.oMi[7].translation()[2] << std::endl;
            }
        }


        // Save trajectory performance for benchmark
        std::ofstream benchFile;
        benchFile.open("analysis/benchmark_data.txt", std::ios_base::app);
        if(benchFile.is_open()){

            // Write extremum of both trajectories
            for(int i=1; i< 29; i++) benchFile << ruckig_traj.row(i).minCoeff() << " ";
            for(int i=1; i< 29; i++) benchFile << ruckig_traj.row(i).maxCoeff() << " ";
            for(int i=1; i< 29; i++) benchFile << polympc_traj.row(i).minCoeff() << " ";
            for(int i=1; i< 29; i++) benchFile << polympc_traj.row(i).maxCoeff() << " ";

            // Write final state to check accuracy
            benchFile << (ruckig_traj.col(nPoints).segment(1,7)-target_position).transpose() << " "
                      << (ruckig_traj.col(nPoints).segment(8,7)-target_velocity).transpose() << " "
                      << (polympc_traj.col(nPoints).segment(1,7)-target_position).transpose() << " "
                      << (polympc_traj.col(nPoints).segment(8,7)-target_velocity).transpose() << " ";
                
            // Write pass/fails flags
            benchFile << jerk_flag_rk << " " 
                      << linear_vel_flag_rk << " " 
                      << angular_vel_flag_rk << " " 
                      << collision_flag_rk << " " 

                      << jerk_flag_mpc << " " 
                      << linear_vel_flag_mpc << " " 
                      << angular_vel_flag_mpc << " " 
                      << collision_flag_mpc << " ";

            // Write target state to check accuracy
            benchFile << (target_position).transpose() << " "
                      << (target_velocity).transpose() << std::endl;
        }

    }

}
