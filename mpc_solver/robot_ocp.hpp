#ifndef SRC_ROCKET_MPC_HPP
#define SRC_ROCKET_MPC_HPP

#include <math.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

#include "polynomials/ebyshev.hpp"
#include "control/continuous_ocp.hpp"
#include "polynomials/splines.hpp"

#include "solvers/sqp_base.hpp"
#include "solvers/osqp_interface.hpp"
#include "control/mpc_wrapper.hpp"

#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/autodiff/cppad.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iostream>
#include <type_traits>
#include <typeinfo>

using CppAD::AD;
using CppAD::NearEqual;

// Global variable with rocket parameters and useful methods

using namespace std;

#define POLY_ORDER 3
#define NUM_SEG    6

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ minTime_ocp, /*NX*/ 14, /*NU*/ 7, /*NP*/ 1, /*ND*/ 0, /*NG*/ 7, /*TYPE*/ double)

using namespace Eigen;

class minTime_ocp : public ContinuousOCP<minTime_ocp, Approximation, SPARSE>{
public:
    
    
    typedef double Scalar;
    typedef AD<Scalar> ADScalar;

    typedef pinocchio::ModelTpl<ADScalar> ADModel;
    typedef ADModel::Data ADData;

    ~minTime_ocp() = default;

    minTime_ocp(){

        pinocchio::urdf::buildModel("robot_utils/panda-model/panda_arm.urdf", model);

        test_ad_model = model.cast<ad_scalar_t>();

        

        // pinocchio::ModelTpl<scalar_t>::Data new_data(model);
        // data = new_data;
    }

    pinocchio::Model model;
    // pinocchio::ModelTpl<scalar_t>::Data data;

    pinocchio::ModelTpl<ad_scalar_t> test_ad_model;

    // ADModel ad_model;

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













   

    EIGEN_STRONG_INLINE void
    inequality_constraints_impl(const Ref<const state_t<scalar_t>> x, const Ref<const control_t<scalar_t>> u,
                                const Ref<const parameter_t <scalar_t>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, Ref<constraint_t <scalar_t>>g) const noexcept
    {

        Eigen::Matrix<double, 7, 1> q = x.head(7);
        Eigen::Matrix<double, 7, 1> q_dot = x.tail(7);
        Eigen::Matrix<double, 7, 1> q_dot_dot = u;

        // Matrix<double, 7, 1> qTarget;

        pinocchio::Data data(model);

        g = pinocchio::rnea(model, data, q, q_dot, q_dot_dot);

    }

    // template<typename T> 
    EIGEN_STRONG_INLINE void
    inequality_constraints_impl(const Ref<const state_t<ad_scalar_t>> x, const Ref<const control_t<ad_scalar_t>> u,
                                const Ref<const parameter_t <ad_scalar_t>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, Ref<constraint_t <ad_scalar_t>>g) const noexcept
    {
        // typedef pinocchio::Model::ConfigVectorType ConfigVectorType;
        // typedef pinocchio::Model::TangentVectorType TangentVectorType;
        // ConfigVectorType q(model.nq);
        // q = pinocchio::randomConfiguration(model);
        // TangentVectorType v(TangentVectorType::Random(model.nv));
        // TangentVectorType a(TangentVectorType::Random(model.nv));

        // typedef ADModel::ConfigVectorType ADConfigVectorType;
        // typedef ADModel::TangentVectorType ADTangentVectorType;

        // ADData ad_data(ad_model);

        // ADConfigVectorType ad_q = q.cast<ADScalar>();
        // ADTangentVectorType ad_v = v.cast<ADScalar>();
        // ADTangentVectorType ad_a = a.cast<ADScalar>();

        // g = pinocchio::rnea(ad_model,ad_data,ad_q,ad_v,ad_a);
    }

    // template<typename T> 
    EIGEN_STRONG_INLINE void
    inequality_constraints_impl(const Ref<const state_t<ad2_scalar_t>> x, const Ref<const control_t<ad2_scalar_t>> u,
                                const Ref<const parameter_t <ad2_scalar_t>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, Ref<constraint_t < ad2_scalar_t>>g) const noexcept
    {
        // g(0) = (T)0.0;
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