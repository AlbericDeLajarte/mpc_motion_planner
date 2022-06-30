#include "pinocchio/codegen/cppadcg.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include <iostream>

// PINOCCHIO_MODEL_DIR is defined by the CMake but you can define your own directory here.
#ifndef PINOCCHIO_MODEL_DIR
  #define PINOCCHIO_MODEL_DIR "path_to_the_model_dir"
#endif

int main(int argc, char ** argv)
{
  using namespace pinocchio;
  
  // You should change here to set up your own URDF file or just pass it as an argument of this example.
//   const std::string urdf_filename = (argc<=1) ? PINOCCHIO_MODEL_DIR + std::string("/example-robot-data/robots/ur_description/urdf/ur5_robot.urdf") : argv[1];

  std::cout << PINOCCHIO_MODEL_DIR << std::endl;
  
  // Load the URDF model
  Model model;
  pinocchio::urdf::buildModel("robot_utils/panda-model/panda_arm.urdf", model);
  // model.gravity = pinocchio::Motion(pinocchio::gravity981());
  model.gravity = pinocchio::Motion::Zero();
  std::cout << model.gravity << std::endl;
  // CodeGenRNEA<double> rnea_code_gen(model);
  // Generate the lib if it does not exist and load it afterwards.
  // rnea_code_gen.initLib();
  // rnea_code_gen.loadLib();

  
  // Build a data related to model
  Data data(model);
  
  // Sample a random joint configuration as well as random joint velocity and acceleration
  Eigen::VectorXd q(7); q << 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0;
  Eigen::VectorXd v = q-q;
  Eigen::VectorXd a(7); a << 10, 10, 10, 10, 10, 10, 10;
  
  // Allocate result container
  Eigen::MatrixXd djoint_torque_dq = Eigen::MatrixXd::Zero(model.nv,model.nv);
  Eigen::MatrixXd djoint_torque_dv = Eigen::MatrixXd::Zero(model.nv,model.nv);
  Eigen::MatrixXd djoint_torque_da = Eigen::MatrixXd::Zero(model.nv,model.nv);

  // Computes the inverse dynamics (RNEA) derivatives for all the joints of the robot
  computeRNEADerivatives(model, data, q, v, a, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);
  crba(model, data, q);
  double lr = 1e-3;

  // std::cout << data.M << std::endl; 
  // std::cout << "-------------\n";

  std::cout << "Initial state:\n" << q.transpose() << std::endl << v.transpose() << std::endl << a.transpose() << std::endl;
  std::cout << "Initial torque norm: " << data.tau.norm() << std::endl;
  std::cout << "theoretical minimum torque: " << rnea(model, data, q, v, q-q).transpose() << std::endl;

  
  for(int i = 0; i< 10000; i++){

    // q -= lr*djoint_torque_dq.transpose()*data.tau;
    // v -= 1*lr*djoint_torque_dv.transpose()*data.tau;
    // a -= 200*lr*djoint_torque_da.transpose()*data.tau;

    computeRNEADerivatives(model, data, q, v, a, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);
    rnea(model, data, q, v, a);
    // data.M.triangularView<Eigen::StrictlyLower>() = data.M.transpose().triangularView<Eigen::StrictlyLower>();

    

    a -= 210*lr*data.M.transpose()*data.tau;

    if (i%1000 == 0) {
      std::cout << "Joint torque: " << data.tau.transpose() << " | Norm = " << data.tau.norm() << std::endl;
      // std::cout << (data.M).transpose() << std::endl;
    }

  }


  std::cout << "Final state:\n" << q.transpose() << std::endl << v.transpose() << std::endl << a.transpose() << std::endl;
  std::cout << "Final torque norm: " << data.tau.norm() << std::endl;
}