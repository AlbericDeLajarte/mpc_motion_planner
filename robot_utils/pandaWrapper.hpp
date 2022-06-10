#pragma once

#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/urdf.hpp"

#define NDOF 7

class PandaWrapper {

    pinocchio::Model model;
    
  public:
    PandaWrapper();
    Eigen::Matrix<double, NDOF, 1> inverse_kinematic(Eigen::Matrix3d orientation, Eigen::Vector3d position);

    std::array<double, NDOF> min_position {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
    std::array<double, NDOF> max_position {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
    std::array<double, NDOF> max_velocity {2.1750, 2.1750, 2.1750, 2.1750, 2.6100, 2.6100, 2.6100};
    std::array<double, NDOF> max_acceleration {15.0, 7.5, 10.0, 12.5, 15.0, 20.0, 20.0};
    std::array<double, NDOF> max_jerk {7500, 3750, 5000, 6250, 7500, 10000, 10000};
};

