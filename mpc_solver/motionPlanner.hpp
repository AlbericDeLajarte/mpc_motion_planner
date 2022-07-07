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
        void warm_start_MPC();

    public:

        MotionPlanner();

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
        void set_current_state(Matrix<double, NDOF, 1> current_position, Matrix<double, NDOF, 1> current_velocity);

        // Set margins on top of robot constraint. Each margin is the ratio of the initial range to be kept
        void set_constraint_margins(double margin_position, double margin_velocity, double margin_acceleration, double margin_torque, double margin_jerk);

        // Return a random position and velocitiy within the bounds (with margin)
        void sample_random_state(Matrix<double, 7, 1> &random_position, Matrix<double, 7, 1> &random_velocity);

        // Solve the OCP to generate the MPC trajectory
        void solve_trajectory();

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

                position_trajectory.col(iPoint) =  Map<Matrix<double, 7, 1> >(new_position.data()),
                velocity_trajectory.col(iPoint) =  Map<Matrix<double, 7, 1> >(new_velocity.data()),
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


        
};