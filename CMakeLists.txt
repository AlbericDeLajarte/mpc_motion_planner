cmake_minimum_required(VERSION 3.0.0)
project(mpc_motion_planner VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pinocchio REQUIRED)
find_package(ruckig REQUIRED)

add_definitions(-DEIGEN_STACK_ALLOCATION_LIMIT=1000000) # 1 MB 

include_directories(submodules/polympc/src mpc_solver robot_utils)
# add_subdirectory(submodules/ruckig)

# add_executable(test_rnea_derivatives test_develop/test_rnea_derivatives.cpp)
# target_link_libraries (test_rnea_derivatives Eigen3::Eigen pinocchio::pinocchio)

add_executable(mpc_benchmark examples/benchmark.cpp robot_utils/pandaWrapper.cpp)
target_link_libraries (mpc_benchmark Eigen3::Eigen pinocchio::pinocchio ruckig::ruckig)

add_executable(offline_trajectory examples/offline_trajectory.cpp robot_utils/pandaWrapper.cpp)
target_link_libraries (offline_trajectory Eigen3::Eigen pinocchio::pinocchio ruckig::ruckig)