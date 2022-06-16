#include "pinocchio/codegen/cppadcg.hpp"
#include "pinocchio/algorithm/crba.hpp"
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

  // CodeGenRNEA<double> rnea_code_gen(model);
  // Generate the lib if it does not exist and load it afterwards.
  // rnea_code_gen.initLib();
  // rnea_code_gen.loadLib();

  
  // Build a data related to model
  Data data(model);
  
  // Sample a random joint configuration as well as random joint velocity and acceleration
  Eigen::VectorXd q(7); q << 1.0, -1.0, 0.5, -2.0, 1.0, 2.0, -0.5;
  Eigen::VectorXd v = 5*q;
  Eigen::VectorXd a = 10*q;
  
  // Allocate result container
  Eigen::MatrixXd djoint_torque_dq = Eigen::MatrixXd::Zero(model.nv,model.nv);
  Eigen::MatrixXd djoint_torque_dv = Eigen::MatrixXd::Zero(model.nv,model.nv);
  Eigen::MatrixXd djoint_torque_da = Eigen::MatrixXd::Zero(model.nv,model.nv);

  // Computes the inverse dynamics (RNEA) derivatives for all the joints of the robot
  computeRNEADerivatives(model, data, q, v, a, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);
  double lr = 0.001;

  std::cout << "Initial state:\n" << q.transpose() << std::endl << v.transpose() << std::endl << a.transpose() << std::endl;
  std::cout << "Initial torque norm: " << data.tau.norm() << std::endl;

  for(int i = 0; i< 100; i++){

    if (i%10 == 0) std::cout << "Joint torque: " << data.tau.transpose() << std::endl;

    // q -= lr*djoint_torque_dq.transpose()*data.tau;
    // v -= 1*lr*djoint_torque_dv.transpose()*data.tau;
    a -= 100*lr*djoint_torque_da.transpose()*data.tau;

    computeRNEADerivatives(model, data, q, v, a, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);

  }

  std::cout << "Final state:\n" << q.transpose() << std::endl << v.transpose() << std::endl << a.transpose() << std::endl;
  std::cout << "Final torque norm: " << data.tau.norm() << std::endl;
    
  // Get access to the joint torque
  
  // std::cout << "torque variation over q:\n " << djoint_torque_dq << std::endl;
  // std::cout << "torque variation over q_dot:\n " << djoint_torque_dv << std::endl;
  // std::cout << "torque variation over q_dot_dot:\n " << djoint_torque_da << std::endl;
}