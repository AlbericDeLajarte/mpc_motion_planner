#include "pandaWrapper.hpp"

PandaWrapper::PandaWrapper(){

    pinocchio::urdf::buildModel("robot_utils/panda-model/panda_arm.urdf", model);

    pinocchio::Data new_data(model);

    data = new_data;
}

Eigen::Matrix<double, 7, 1> PandaWrapper::inverse_kinematic(Eigen::Matrix3d orientation, Eigen::Vector3d position){

    const pinocchio::SE3 oMdes(orientation, position);
    
    Eigen::VectorXd q = pinocchio::randomConfiguration(model);
    const double eps  = 1e-4;
    const int IT_MAX  = 1000;
    const double DT   = 1e-1;
    const double damp = 1e-2;
    
    pinocchio::Data::Matrix6x J(6,model.nv);
    J.setZero();
    
    bool success = false;
    Eigen::Matrix<double, 6, 1> err;
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
    }

    Eigen::Ref<Eigen::Matrix<double, 7, 1>> qSolution(q);

    return qSolution;

}

Eigen::Matrix<double, NDOF, 1> PandaWrapper::inverse_velocities(Eigen::Matrix<double, NDOF, 1> q, Eigen::Vector3d linear_velocity, Eigen::Vector3d angular_velocity){

    // Construct data
    pinocchio::Data::Matrix6x J(6, model.nv);
    J.setZero();
    pinocchio::computeJointJacobian(model, data, q, JOINT_ID, J);

    Eigen::VectorXd joint_velocity(model.nv);
    Eigen::Matrix<double, 6, 1> task_velocity;
    task_velocity << linear_velocity, angular_velocity; 
    
    // Solve using damped pseudo-inverse
    pinocchio::Data::Matrix6 JJt;
    JJt.noalias() = J * J.transpose();
    JJt.diagonal().array() += 1e-3;
    joint_velocity.noalias() = J.transpose() * JJt.ldlt().solve(task_velocity);

    // Cast to static matrix
    Eigen::Ref<Eigen::Matrix<double, NDOF, 1>> joint_velocity_solution(joint_velocity);
    return joint_velocity_solution;
}