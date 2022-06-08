cmake_minimum_required(VERSION 3.0.0)
project(mpc_motion_planner VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(mpc_motion_planner main.cpp)
target_link_libraries (mpc_motion_planner Eigen3::Eigen)

include_directories(submodules/polympc/src)