#include "main.hpp"
#include "pinocchio/algorithm/rnea.hpp"



int main(int, char**) { 

    std::srand((unsigned int) time(0));

    MotionPlanner planner;

    // Add margins on limits
    planner.robot.min_position = planner.robot.min_position.array()+0.3;
    planner.robot.max_position = planner.robot.max_position.array()-0.3;

    planner.robot.max_velocity *= 0.5;
    // planner.robot.max_velocity  = planner.robot.max_velocity.array() -0.3;
    planner.robot.max_acceleration *= 0.1;
    planner.robot.max_jerk *= 0.01;
    planner.robot.max_torque *= 0.7;
    // planner.robot.max_torqueDot *= 0.7;
    

    // ---------- Solve trajectories ---------- //
    
    for(int iter =0; iter<2; iter++){
        
        
        Matrix<double, 7, 1> target_position = 0.5*(Matrix<double, 7, 1>::Random().array()*(planner.robot.max_position - planner.robot.min_position).array() 
                                                    + (planner.robot.max_position + planner.robot.min_position).array() );
        // Matrix<double, 7, 1> target_velocity = Matrix<double, 7, 1>::Random().array()*plamner.robot.max_velocity.array();
        Matrix<double, 7, 1> target_velocity = planner.robot.inverse_velocities(target_position, Matrix<double, 3, 1>::Random()*planner.robot.max_linear_velocity, Vector3d(0.0, 0.0, 0.0));

        Matrix<double, 6, 1> task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
        
        // std::cout << "position: " << target_position.transpose() << std::endl;
        // std::cout << "velocity: " << target_velocity.transpose() << std::endl;
        // std::cout << "jacobian:\n" << J << std::endl;
        // std::cout << "task velocity: " << task_velocity.transpose() << std::endl;
        
        if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) {
            target_velocity *= 0.9 * planner.robot.max_linear_velocity / task_velocity.head(3).norm();
            
            std::cout << "Linear Vel: " << task_velocity.head(3).norm();
            task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
            std::cout << " corrected to : " << task_velocity.head(3).norm()  << std::endl;
        }
        if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) {
            target_velocity *= 0.9 * planner.robot.max_angular_velocity / task_velocity.tail(3).norm();
            
            std::cout << "Angular Vel: " << task_velocity.tail(3).norm();
            task_velocity = planner.robot.forward_velocities(target_position, target_velocity);
            std::cout << " corrected to : " << task_velocity.tail(3).norm() << std::endl;
        }

        if( (target_velocity.array().abs() > planner.robot.max_velocity.array()).any() ){ 
            // std::cout << "RK: Jerk limit of: " << jerk.transpose() << " at time: " << time << std::endl;
            target_velocity /= 1.1*(target_velocity.array().abs()/planner.robot.max_velocity.array()).maxCoeff();
        }

                
        planner.set_target_state(target_position, target_velocity);
        std::cout << planner.current_state.reshaped(7, 2).transpose() << std::endl; 
        std::cout << planner.target_state.reshaped(7, 2).transpose() << std::endl; 


        
        // ---------- SOLVE POLYMPC ---------- //
        planner.solve_trajectory();

        

        // --------- Write data to txt file --------- //
        const int nPoints = 200;

        int linear_vel_flag_rk {1}, angular_vel_flag_rk {1}, jerk_flag_rk {1}, torqueDot_flag_rk {1};
        int linear_vel_flag_mpc {1}, angular_vel_flag_mpc {1}, jerk_flag_mpc {1}, torqueDot_flag_mpc {1};


        // Check Ruckig trajectory
        Matrix<double, 7, nPoints+1> position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory;
        planner.get_ruckig_trajectory<nPoints>(position_trajectory, velocity_trajectory, acceleration_trajectory, torque_trajectory);
        
        double dT = planner.trajectory.get_duration()/nPoints;
        for (int iPoint = 0; iPoint<=nPoints; iPoint++)
        {
            double time = dT * iPoint;

            // Check jerk variation
            if (iPoint >= 1){
                Matrix<double, 7, 1> jerk = (acceleration_trajectory.col(iPoint)-acceleration_trajectory.col(iPoint-1)) / dT;

                if( (jerk.array().abs() > 10*planner.robot.max_jerk.array()).any() ){ // Using 10% margin because RK is driving at the limit
                    jerk_flag_rk = 0;
                    // std::cout << "RK: Jerk limit of: " << jerk.transpose() << " at time: " << time << std::endl;
                } 
            }

            // Check cartesian velocity
            task_velocity = planner.robot.forward_velocities(position_trajectory.col(iPoint), velocity_trajectory.col(iPoint));

            if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) {
                linear_vel_flag_rk = 0;
                std::cout << "RK: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << time << std::endl;
            }
            if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) {
                angular_vel_flag_rk = 0;
                std::cout << "RK: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << time << std::endl;
            }
        }
        

            // // Log Polympc trajectory
            // dT = mpc.solution_p()[0]/nPoints;
            // for (int iPoint = 0; iPoint<=nPoints; iPoint++)
            // {
            //     double time = 1.0/nPoints * iPoint;

            //     polympc_traj.col(iPoint).head(21) << mpc.solution_x_at(time), mpc.solution_u_at(time);

            //     polympc_traj.col(iPoint).tail(7) = 
            //     pinocchio::rnea(planner.robot.model, planner.robot.data, polympc_traj.col(iPoint).head(7),
            //                                             polympc_traj.col(iPoint).segment(7, 7),
            //                                             polympc_traj.col(iPoint).segment(14, 7));
                
            //     // Check jerk and torque variation
            //     if (iPoint >= 1){
            //         Matrix<double, 7, 1> jerk = (polympc_traj.col(iPoint).segment(14, 7)-polympc_traj.col(iPoint-1).segment(14, 7)) / dT;
            //         Matrix<double, 7, 1> torqueDot = (polympc_traj.col(iPoint).tail(7)-polympc_traj.col(iPoint-1).tail(7)) / dT;

            //         if( (jerk.array().abs() > 10*planner.robot.max_jerk.array()).any() ){
            //             jerk_flag_mpc = 0;
            //             // std::cout << "MPC: Jerk limit of: " << jerk.transpose() << " at time: " << time << std::endl;
            //         } 
            //         if( (torqueDot.array().abs() > planner.robot.max_torqueDot).any() ){
            //             torqueDot_flag_mpc = 0;
            //             // std::cout << "MPC: TorqueDot limit of: " << torqueDot.transpose() << " at time: " << time << std::endl;
            //         } 
            //     }

            //     // Check cartesian velocity
            //     task_velocity = planner.robot.forward_velocities(polympc_traj.col(iPoint).head(7), polympc_traj.col(iPoint).segment(7, 7));

            //     if(task_velocity.head(3).norm() > planner.robot.max_linear_velocity) {
            //         linear_vel_flag_mpc = 0;
            //         std::cout << "MPC: Linear Vel limit of: " << task_velocity.head(3).norm() << " at time: " << time << std::endl;
            //     }
            //     if(task_velocity.tail(3).norm() > planner.robot.max_angular_velocity) {
            //         angular_vel_flag_mpc = 0;
            //         std::cout << "MPC: Angular Vel limit of: " << task_velocity.tail(3).norm() << " at time: " << time << std::endl;
            //     }

            //     logFile << time*mpc.solution_p()[0] << " " 
            //             << polympc_traj.col(iPoint).transpose() 
            //             << std::endl;
            // }


        // Save trajectory performance for benchmark
        // std::ofstream benchFile;
        // benchFile.open("data/benchmark_data.txt", std::ios_base::app);
        // if(benchFile.is_open()){

        //     // Log Ruckig data  
        //     trajectory.at_time(trajectory.get_duration(), new_position, new_velocity, new_acceleration);

        //     // Write extremum of both trajectories
        //     for(int i=0; i< 28; i++) benchFile << ruckig_traj.row(i).minCoeff() << " ";
        //     for(int i=0; i< 28; i++) benchFile << ruckig_traj.row(i).maxCoeff() << " ";
        //     for(int i=0; i< 28; i++) benchFile << polympc_traj.row(i).minCoeff() << " ";
        //     for(int i=0; i< 28; i++) benchFile << polympc_traj.row(i).maxCoeff() << " ";

        //     // Write final state
        //     benchFile << (ruckig_traj.col(nPoints).head(14)-final_state).transpose() << " "
        //             << (polympc_traj.col(nPoints).head(14)-final_state).transpose() << " ";
                
        //     benchFile << jerk_flag_rk << " " 
        //             << torqueDot_flag_rk << " " 
        //             << linear_vel_flag_rk << " " 
        //             << angular_vel_flag_rk << " " 

        //             << jerk_flag_mpc << " " 
        //             << torqueDot_flag_mpc << " "
        //             << linear_vel_flag_mpc << " " 
        //             << angular_vel_flag_mpc << std::endl;
        // }

    }

}
