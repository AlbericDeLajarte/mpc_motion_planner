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

    public:

        MotionPlanner(): robot(), mpc(){

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

            // State constraints ---------------
            mpc_t::state_t lbx; lbx << robot.min_position, -robot.max_velocity;
            mpc_t::state_t ubx; ubx << robot.max_position, robot.max_velocity;
            mpc.state_bounds(lbx, ubx);

            // Input constraints -------------
            mpc.control_bounds(-robot.max_acceleration, robot.max_acceleration);  
            
            // Parameters ------------------
            mpc_t::parameter_t lbp; lbp << 0.0;  // lower bound on time
            mpc_t::parameter_t ubp; ubp << 10;   // upper bound on time
            mpc.parameters_bounds(lbp, ubp);

            // Non-linear torque constraints
            mpc.constraints_bounds(-robot.max_torque, robot.max_torque);

            // ---------- Ruckig setup ---------- //
            Matrix<double, 7, 1>::Map(input.max_velocity.data() ) = robot.max_velocity;
            Matrix<double, 7, 1>::Map(input.max_acceleration.data() ) = robot.max_acceleration;
            Matrix<double, 7, 1>::Map(input.max_jerk.data() ) = robot.max_jerk;
            
        }

        // Creates solver
        using mpc_t = MPC<minTime_ocp, MySolver, admm>;
        mpc_t mpc;

        PandaWrapper robot;

        // Ruckig as warm start
        Ruckig<NDOF> otg;
        Trajectory<NDOF> trajectory;
        InputParameter<NDOF> input;

        // Robot state and target
        Eigen::Matrix<double, NDOF, 1> init_position; // to delete
        Eigen::Matrix<double, NDOF, 1> init_velocity; // to delete

        Eigen::Matrix<double, 2*NDOF, 1> current_state;
        Eigen::Matrix<double, 2*NDOF, 1> target_state;

        const double eps = 1e-2;
        const double inf = std::numeric_limits<double>::infinity();
        

        void set_target_state(Matrix<double, NDOF, 1> target_position, Matrix<double, NDOF, 1> target_velocity){

            // Update MPC constraints
            target_state.head(7) = target_position;
            target_state.tail(7) = target_velocity;

            mpc.final_state_bounds(target_state.array() - eps, target_state.array() + eps);

            // Update Ruckig constraints
            Matrix<double, 7, 1>::Map(input.target_position.data() ) = target_position;
            Matrix<double, 7, 1>::Map(input.target_velocity.data() ) = target_velocity;
            input.target_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  

            // Warm start when we change target = new trajectory
            warm_start_MPC();
        }

        void set_current_state(Matrix<double, NDOF, 1> current_position, Matrix<double, NDOF, 1> current_velocity){

            // Update MPC constraints
            current_state.head(7) = current_position;
            current_state.tail(7) = current_velocity;

            mpc.initial_conditions(current_state);


            // Update Ruckig constraints
            Matrix<double, 7, 1>::Map(input.current_position.data() ) = current_position;
            Matrix<double, 7, 1>::Map(input.current_velocity.data() ) = current_velocity;
            input.current_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};    
        }

        void warm_start_MPC(){

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

        void solve_trajectory(){
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