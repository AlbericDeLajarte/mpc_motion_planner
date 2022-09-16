#ifndef SRC_ROCKET_MPC_HPP
#define SRC_ROCKET_MPC_HPP

#include <math.h>

#include <Eigen/Dense>

#include "polynomials/ebyshev.hpp"
#include "control/continuous_ocp.hpp"
#include "polynomials/splines.hpp"

#include "solvers/sqp_base.hpp"
#include "solvers/osqp_interface.hpp"
#include "control/mpc_wrapper.hpp"

#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
// #include "robotDynamic.hpp"
// #include "pinocchio/algorithm/joint-configuration.hpp"

#include <iostream>

// Global variable with rocket parameters and useful methods

using namespace std;

#define POLY_ORDER 3
#define NUM_SEG    6

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ minTime_ocp, /*NX*/ 14, /*NU*/ 7, /*NP*/ 1, /*ND*/ 0, /*NG*/ 8, /*TYPE*/ double)

using namespace Eigen;

class minTime_ocp : public ContinuousOCP<minTime_ocp, Approximation, SPARSE>{
public:

    ~minTime_ocp() = default;

    pinocchio::Model model;

    void init(std::string urdf_path){
        pinocchio::urdf::buildModel(urdf_path, model);
    }

    template<typename T>
    inline void dynamics_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                              const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> &d,
                              const T &t, Eigen::Ref<state_t<T>> xdot) const noexcept
    {

        // -------------- Differential equation ---------------------

        // Position variation is joint velocity
        xdot.head(7) = x.segment(7, 7);

        // Joint velocity variation is input = acceleration
        xdot.segment(7, 7) = u;

        // Scaling dynamic with time parameter        
        xdot *= p(0);
        
        polympc::ignore_unused_var(t);
    }






EIGEN_STRONG_INLINE constraint_t<scalar_t> evalConstraints(const Ref<const state_t<scalar_t>> x, const Ref<const control_t<scalar_t>> u) const noexcept
{
    Eigen::Matrix<double, 7, 1> q = x.head(7);
    Eigen::Matrix<double, 7, 1> q_dot = x.tail(7);
    Eigen::Matrix<double, 7, 1> q_dot_dot = u;

    pinocchio::Data data(model);
    pinocchio::forwardKinematics(model, data, q);

    constraint_t<scalar_t> ineq_constraint;
    ineq_constraint << pinocchio::rnea(model, data, q, q_dot, q_dot_dot),  data.oMi[7].translation()[2];

    // std::cout << "----\n" << ineq_constraint.transpose() << "\n----\n";

    return ineq_constraint;
}

EIGEN_STRONG_INLINE constraint_t<ad_scalar_t> evalConstraints(const Ref<const state_t<ad_scalar_t>> x, const Ref<const control_t<ad_scalar_t>> u) const noexcept
{
    Eigen::Matrix<double, 7, 1> q;
    Eigen::Matrix<double, 7, 1> q_dot;
    Eigen::Matrix<double, 7, 1> q_dot_dot;

    for(int i = 0; i<7; i++){
        q(i) = x(i).value();
        q_dot(i) = x(i+7).value();
        q_dot_dot(i) = u(i).value();
    }
    

    pinocchio::Data data(model);

    // Allocate result container
    Eigen::MatrixXd djoint_torque_dq = Eigen::MatrixXd::Zero(model.nv,model.nv);
    Eigen::MatrixXd djoint_torque_dv = Eigen::MatrixXd::Zero(model.nv,model.nv);
    Eigen::MatrixXd djoint_torque_da = Eigen::MatrixXd::Zero(model.nv,model.nv);

    pinocchio::computeRNEADerivatives(model, data, q, q_dot, q_dot_dot, djoint_torque_dq, djoint_torque_dv, djoint_torque_da);

    pinocchio::rnea(model, data, q, q_dot, q_dot_dot);
    pinocchio::crba(model, data, q);
    data.M.triangularView<Eigen::StrictlyLower>() = data.M.transpose().triangularView<Eigen::StrictlyLower>();

    Eigen::MatrixXd djoint_torque_dtime_f = djoint_torque_dv*q_dot + djoint_torque_da*q_dot_dot;

    constraint_t<ad_scalar_t> ineq_constraint;

    // Fill in torque derivatives
    Eigen::Matrix<scalar_t, 1, NX + NU + NP> jac_row;
    jac_row.setZero();
    for(int i = 0; i < NG-1; i++)
    {   
        // Overwitting with PINOCCHIO data
        jac_row.head(7) = djoint_torque_dq.row(i);
        jac_row.segment(7, 7) = djoint_torque_dv.row(i);
        jac_row.segment(14, 7) = data.M.row(i);

        jac_row(21) = djoint_torque_dtime_f(i);

        ineq_constraint(i).value() = data.tau(i);
        ineq_constraint(i).derivatives() = jac_row;
    }

    // Fill in height derivative
    pinocchio::Data::Matrix6x J(6,7); J.setZero();

    pinocchio::forwardKinematics(model, data, q);
    pinocchio::computeJointJacobian(model, data, q, 7, J);

    // For some reasons we should rotate the jacobian
    Eigen::MatrixXd rotateJacobian(6,6); rotateJacobian.setZero();
    rotateJacobian.block(0,0,3,3) = data.oMi[7].rotation();
    rotateJacobian.block(3,3,3,3) = data.oMi[7].rotation();
    J = rotateJacobian * J;

    ineq_constraint(NG-1).value() = data.oMi[7].translation()[2];
    jac_row = Eigen::Matrix<scalar_t, 1, NX + NU + NP>::Zero();
    jac_row.head(7) = J.row(2);
    ineq_constraint(NG-1).derivatives() = jac_row;

    return ineq_constraint;
}

EIGEN_STRONG_INLINE constraint_t<ad2_scalar_t> evalConstraints(const Ref<const state_t<ad2_scalar_t>> x, const Ref<const control_t<ad2_scalar_t>> u) const noexcept
{
    return constraint_t<ad2_scalar_t>::Zero();
}


template<typename T>
EIGEN_STRONG_INLINE void
inequality_constraints_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
                            const Ref<const parameter_t <T>> p, const Ref<const static_parameter_t> d,
                            const scalar_t &t, Ref<constraint_t < T>>g) const noexcept
{
    g = evalConstraints(x, u);
}
    













    template<typename T>
    inline void lagrange_term_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                                   const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> d,
                                   const scalar_t &t, T &lagrange) noexcept
    {
        lagrange = (T)0.0;
    }

    template<typename T>
    inline void mayer_term_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                                const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> d,
                                const scalar_t &t, T &mayer) noexcept
    {   
        mayer = p(0);

        // polympc::ignore_unused_var(x);
        polympc::ignore_unused_var(u);
        polympc::ignore_unused_var(d);
        polympc::ignore_unused_var(t);
  
    }
};

#endif //SRC_ROCKET_MPC_HPP