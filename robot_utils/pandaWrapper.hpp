#pragma once

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/urdf.hpp"

#define NDOF 7

class PandaWrapper {

    
    
  public:

    pinocchio::Model model;
    pinocchio::Data data;
    int frame_id;

    PandaWrapper(std::string urdf_path);
    Eigen::Matrix<double, NDOF, 1> inverse_kinematic(Eigen::Matrix3d orientation, Eigen::Vector3d position);
    Eigen::Matrix<double, NDOF, 1> inverse_velocities(Eigen::Matrix<double, NDOF, 1> q, Eigen::Vector3d linear_velocity, Eigen::Vector3d angular_velocity);
    Eigen::Matrix<double, 6, 1> forward_velocities(Eigen::Matrix<double, NDOF, 1> q, Eigen::Matrix<double, NDOF, 1> qdot);


    // Limits from https://frankaemika.github.io/docs/control_parameters.html
    Eigen::Matrix<double, NDOF, 1> min_position {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
    Eigen::Matrix<double, NDOF, 1> max_position {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
    Eigen::Matrix<double, NDOF, 1> max_velocity {2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100};
    Eigen::Matrix<double, NDOF, 1> max_acceleration {15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0};
    Eigen::Matrix<double, NDOF, 1> max_jerk {7500, 3750, 5000, 6250, 7500, 10000, 10000};
    Eigen::Matrix<double, NDOF, 1> max_torque {87, 87, 87, 87, 12, 12, 12};
    double max_torqueDot {1000};

    double max_linear_velocity {1.7};
    double max_angular_velocity {2.5};

    double min_height {0.05};
    
};

