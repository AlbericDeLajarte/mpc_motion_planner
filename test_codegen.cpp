//
 // Copyright (c) 2020 INRIA
 //
  
 #include "pinocchio/codegen/cppadcg.hpp" // this file should be included first before all the others!
 #include "pinocchio/algorithm/rnea.hpp"
  
 #include "pinocchio/parsers/urdf.hpp"
 #include "pinocchio/algorithm/joint-configuration.hpp"
 #include "pinocchio/codegen/code-generator-algo.hpp"
  
 #include <iostream>
  
 int main(int argc, const char ** argv)
 {
   using namespace pinocchio;
   using namespace Eigen;
   
      
   // Load the model
   Model model; pinocchio::urdf::buildModel("robot_utils/panda-model/panda_arm.urdf", model);
    pinocchio::Data data(model);

   CodeGenRNEA<double> rnea_code_gen(model, "rnea", "cg_rnea");
  
   // Generate the lib if it does not exist and load it afterwards.
   rnea_code_gen.initLib();
   rnea_code_gen.loadLib();
  
   // Use it with a random configuration samples in the bounds of the joint limits
   // Sample a random joint configuration as well as random joint velocity and acceleration

   for (int i = 0; i<20; i++){
        VectorXd q = randomConfiguration(model);
        VectorXd v = 5*q;
        VectorXd a = 10*q;
        rnea_code_gen.evalFunction(q, v, a);


        // Retrieve the result
        std::cout << rnea_code_gen.rnea(q, v, a).transpose() << std::endl;
        std::cout << rnea(model, data, q, v, a).transpose() << std::endl;
   }
  
   // And make it symmetric if needed
//    M.template triangularView<Eigen::StrictlyLower>() = M.transpose().template triangularView<Eigen::StrictlyLower>();
  
//    // You can check the result with the classic rnea
//    Data data_check(model);
//    rnea(model,data_check,q);
  
//    data_check.M.triangularView<Eigen::StrictlyLower>() = data_check.M.transpose().triangularView<Eigen::StrictlyLower>();
  
//    const MatrixXd & M_check = data_check.M;
//    if(M_check.isApprox(M)) {
//      std::cout << "Super! The two results are the same." << std::endl;
//      return 0;
//    }
//    else {
//      std::cout << "Not Super! The results do not match." << std::endl;
//      return -1;
//    }
  
 }