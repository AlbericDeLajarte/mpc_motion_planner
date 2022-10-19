#include <iostream>
#include <fstream>
#include "robot_ocp.hpp"
#include "polympc_redef.hpp" 

#include "pandaWrapper.hpp"

#include <ruckig/ruckig.hpp>

using admm = boxADMM<minTime_ocp::VAR_SIZE, minTime_ocp::NUM_EQ + minTime_ocp::NUM_INEQ, minTime_ocp::scalar_t,
                minTime_ocp::MATRIXFMT, linear_solver_traits<minTime_ocp::MATRIXFMT>::default_solver>;

using namespace Eigen;
using namespace ruckig;

class MotionPlanner{

    private:

        // Solve problem with ruckig to initialize MPC
        void warm_start_RK();

    public:

        MotionPlanner(std::string urdf_path);

        // PolyMPC solver
        using mpc_t = MPC<minTime_ocp, MySolver, admm>;
        mpc_t mpc;

        // Utility class wrapping pinocchio
        PandaWrapper robot;

        // Ruckig as warm start
        Ruckig<NDOF> otg;
        Trajectory<NDOF> trajectory;
        InputParameter<NDOF> input;

        // Robot state and target
        Eigen::Matrix<double, 2*NDOF, 1> current_state;
        Eigen::Matrix<double, 2*NDOF, 1> target_state;

        // Utility attributes for solver constraints
        const double eps = 1e-2;
        const double inf = std::numeric_limits<double>::infinity();

        // -------- Methods -------- //

        // Margins on bounds
        double margin_position_, margin_velocity_, margin_acceleration_, margin_torque_, margin_jerk_;
        
        // Set the target (final state) as a constraint
        void set_target_state(Matrix<double, NDOF, 1> target_position, Matrix<double, NDOF, 1> target_velocity);

        // Set the current (initial state) as a constraint
        void set_current_state(Matrix<double, NDOF, 1> current_position, Matrix<double, NDOF, 1> current_velocity, Matrix<double, NDOF, 1> current_acceleration = Matrix<double, NDOF, 1>::Zero());

        // Set margins on top of robot constraint. Each margin is the ratio of the initial range to be kept
        void set_constraint_margins(double margin_position, double margin_velocity, double margin_acceleration, double margin_torque, double margin_jerk);

        void set_min_height(double min_height);
        
        // Return a random position and velocitiy within the bounds (with margin)
        void sample_random_state(Matrix<double, 7, 1> &random_position, Matrix<double, 7, 1> &random_velocity);

        // Solve the OCP to generate the MPC trajectory
        void solve_trajectory(bool use_ruckig_as_warm_start);

        // Check if a given state is feasible
        int check_state_in_bounds(Matrix<double, 7, 1> &position, Matrix<double, 7, 1> &velocity);

        // Get the N points from ruckig trajectory
        template<const int N>
        void get_ruckig_trajectory(Matrix<double, 1, N+1> &time, Matrix<double, 7, N+1> &position_trajectory, Matrix<double, 7, N+1> &velocity_trajectory, Matrix<double, 7, N+1> &acceleration_trajectory, Matrix<double, 7, N+1> &torque_trajectory){

            // Check Ruckig trajectory
            double dT = trajectory.get_duration()/N;

            std::array<double, NDOF> new_position, new_velocity, new_acceleration;
            
            for (int iPoint = 0; iPoint<=N; iPoint++)
            {
                time(iPoint) = dT * iPoint;
                
                trajectory.at_time(time(iPoint), new_position, new_velocity, new_acceleration);

                position_trajectory.col(iPoint) =  Map<Matrix<double, 7, 1> >(new_position.data());
                velocity_trajectory.col(iPoint) =  Map<Matrix<double, 7, 1> >(new_velocity.data());
                acceleration_trajectory.col(iPoint) = Map<Matrix<double, 7, 1> >(new_acceleration.data());

                torque_trajectory.col(iPoint) =
                pinocchio::rnea(robot.model, robot.data, position_trajectory.col(iPoint),
                                                            velocity_trajectory.col(iPoint),
                                                            acceleration_trajectory.col(iPoint));
            }
        }

        // Get the N points from MPC trajectory
        template<const int N>
        void get_MPC_trajectory(Matrix<double, 1, N+1> &time, Matrix<double, 7, N+1> &position_trajectory, Matrix<double, 7, N+1> &velocity_trajectory, Matrix<double, 7, N+1> &acceleration_trajectory, Matrix<double, 7, N+1> &torque_trajectory){

            for (int iPoint = 0; iPoint<=N; iPoint++)
            {
                time(iPoint) = 1.0/N * iPoint;
                    
                position_trajectory.col(iPoint) = mpc.solution_x_at(time(iPoint)).head(7);
                velocity_trajectory.col(iPoint) = mpc.solution_x_at(time(iPoint)).tail(7);
                acceleration_trajectory.col(iPoint) = mpc.solution_u_at(time(iPoint));

                torque_trajectory.col(iPoint) = 
                pinocchio::rnea(robot.model, robot.data, position_trajectory.col(iPoint),
                                                            velocity_trajectory.col(iPoint),
                                                            acceleration_trajectory.col(iPoint));
            }
            time *= mpc.solution_p()[0];
        }

        void get_MPC_point(double time, Matrix<double, 7, 1> &position, Matrix<double, 7, 1> &velocity, Matrix<double, 7, 1> &acceleration, Matrix<double, 7, 1> &torque){
            
            if (time < mpc.solution_p()[0]) time /= mpc.solution_p()[0];
            else time = mpc.solution_p()[0];

            position = mpc.solution_x_at(time).head(7);
            velocity = mpc.solution_x_at(time).tail(7);
            acceleration = mpc.solution_u_at(time);

            torque = pinocchio::rnea(robot.model, robot.data, position, velocity, acceleration);
        }

        void get_RK_point(double time, Matrix<double, 7, 1> &position, Matrix<double, 7, 1> &velocity, Matrix<double, 7, 1> &acceleration, Matrix<double, 7, 1> &torque){

            time = std::min(time, trajectory.get_duration());

            std::array<double, NDOF> new_position, new_velocity, new_acceleration;
            trajectory.at_time(time, new_position, new_velocity, new_acceleration);

            position = Map<Matrix<double, 7, 1> >(new_position.data());
            velocity = Map<Matrix<double, 7, 1> >(new_velocity.data());
            acceleration = Map<Matrix<double, 7, 1> >(new_acceleration.data());

            torque = pinocchio::rnea(robot.model, robot.data, position, velocity, acceleration);
        }

        // Warm start MPC with given trajectory (should be regularly time spaced)
        void warm_start(double final_time, MatrixXd position_trajectory, MatrixXd velocity_trajectory, MatrixXd acceleration_trajectory){
            
            int nPoint = position_trajectory.cols();

            mpc_t::traj_state_t x_guess;
            mpc_t::traj_control_t u_guess;
            mpc_t::parameter_t p0; p0 << final_time;

            auto mpc_time_grid = mpc.ocp().time_nodes;
            int i = 0;
            for(auto mpc_time : mpc_time_grid){

                int traj_idx = std::round(mpc_time*(nPoint-1));

                x_guess.segment(i*NDOF*2, NDOF*2) << position_trajectory.col(traj_idx),
                                                     velocity_trajectory.col(traj_idx);

                u_guess.segment(i*NDOF, NDOF) << acceleration_trajectory.col(traj_idx);;

                i++;
            } 
            // std::cout << x_guess << std::endl << x_guess.cols() << " " << x_guess.rows() << std::endl;
            // std::cout << x_guess.reshaped(14, 13)  << std::endl;
            // std::cout << u_guess.reshaped(7, 13)  << std::endl;
            mpc.x_guess(x_guess);	
            mpc.u_guess(u_guess);
            mpc.p_guess(p0); 
        }


        
};