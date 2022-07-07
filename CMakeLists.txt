cmake_minimum_required(VERSION 3.0.0)
project(mpc_motion_planner_examples VERSION 0.1.0)

set(CMAKE_BUILD_TYPE Release)

add_subdirectory(mpc_solver)

# add_executable(test_rnea_derivatives test_develop/test_rnea_derivatives.cpp)
# target_link_libraries (test_rnea_derivatives Eigen3::Eigen pinocchio::pinocchio)

add_executable(mpc_benchmark examples/benchmark.cpp)
target_link_libraries (mpc_benchmark PUBLIC motionPlanner)

add_executable(offline_trajectory examples/offline_trajectory.cpp)
target_link_libraries (offline_trajectory PUBLIC motionPlanner)