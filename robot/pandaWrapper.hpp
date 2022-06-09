#pragma once

#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/urdf.hpp"

class PandaWrapper {

    pinocchio::Model model;
    
  public:
    PandaWrapper();
    Eigen::VectorXd inverse_kinematic(Eigen::Matrix3d orientation, Eigen::Vector3d position);
};

