#include "pandaWrapper.hpp"

PandaWrapper::PandaWrapper(){

    pinocchio::urdf::buildModel("panda-model/panda_arm.urdf", model);
        
}

Eigen::VectorXd PandaWrapper::inverse_kinematic(Eigen::Matrix3d orientation, Eigen::Vector3d position){

    pinocchio::Data data(model);

    const int JOINT_ID = 7;
    const pinocchio::SE3 oMdes(orientation, position);
    
    Eigen::VectorXd q = pinocchio::neutral(model);
    const double eps  = 1e-4;
    const int IT_MAX  = 1000;
    const double DT   = 1e-1;
    const double damp = 1e-6;
    
    pinocchio::Data::Matrix6x J(6,model.nv);
    J.setZero();
    
    bool success = false;
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    Vector6d err;
    Eigen::VectorXd v(model.nv);
    for (int i=0;;i++)
    {
        pinocchio::forwardKinematics(model,data,q);
        const pinocchio::SE3 dMi = oMdes.actInv(data.oMi[JOINT_ID]);
        err = pinocchio::log6(dMi).toVector();
        if(err.norm() < eps)
        {
        success = true;
        break;
        }
        if (i >= IT_MAX)
        {
        success = false;
        break;
        }
        pinocchio::computeJointJacobian(model,data,q,JOINT_ID,J);
        pinocchio::Data::Matrix6 JJt;
        JJt.noalias() = J * J.transpose();
        JJt.diagonal().array() += damp;
        v.noalias() = - J.transpose() * JJt.ldlt().solve(err);
        q = pinocchio::integrate(model,q,v*DT);
        if(!(i%10))
        std::cout << i << ": error = " << err.transpose() << std::endl;
    }
    
    if(success) 
    {
        std::cout << "Convergence achieved!" << std::endl;
    }
    else 
    {
        std::cout << "\nWarning: the iterative algorithm has not reached convergence to the desired precision" << std::endl;
    }
        
    std::cout << "\nresult: " << q.transpose() << std::endl;
    std::cout << "\nfinal error: " << err.transpose() << std::endl;

    return q;

}