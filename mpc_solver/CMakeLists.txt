cmake_minimum_required(VERSION 3.0.0)
project(mpc_motion_planner VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pinocchio REQUIRED)
find_package(ruckig REQUIRED)

add_definitions(-DEIGEN_STACK_ALLOCATION_LIMIT=1000000) # 1 MB 

include_directories(../submodules/polympc/src ../robot_utils)

add_library(motionPlanner motionPlanner.cpp ../robot_utils/pandaWrapper.cpp)
target_link_libraries (motionPlanner Eigen3::Eigen pinocchio::pinocchio ruckig::ruckig)
target_include_directories(motionPlanner PUBLIC "." "../submodules/polympc/src" "../robot_utils")