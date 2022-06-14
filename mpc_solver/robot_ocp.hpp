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

#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include <iostream>
#include <type_traits>
#include <typeinfo>


// Global variable with rocket parameters and useful methods

using namespace std;

#define POLY_ORDER 3
#define NUM_SEG    6

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ minTime_ocp, /*NX*/ 14, /*NU*/ 7, /*NP*/ 1, /*ND*/ 0, /*NG*/ 1, /*TYPE*/ double)

using namespace Eigen;

class minTime_ocp : public ContinuousOCP<minTime_ocp, Approximation, SPARSE>{
public:
    ~minTime_ocp() = default;

    minTime_ocp(){

        pinocchio::urdf::buildModel("robot_utils/panda-model/panda_arm.urdf", model);

        pinocchio::Data new_data(model);

        data = new_data;
    }

    pinocchio::Model model;
    pinocchio::Data data;

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













    // template<typename T>
    // typename std::enable_if<T == scalar_t>::type
    // template < typename = typename std::enable_if< true >::type >
    // inline void inequality_constraints_impl(const Ref<const state_t<scalar_t>> x, const Ref<const control_t<scalar_t>> u,
    //                                         const Ref<const parameter_t<scalar_t>> p, const Ref<const static_parameter_t> d,
    //                                         const scalar_t &t, Eigen::Ref<constraint_t<scalar_t>> g) const noexcept
    // {
    //     Matrix<double, 7, 1> q;
    //     Matrix<double, 7, 1> q_dot;
    //     Matrix<double, 7, 1> q_dot_dot;
    //     g(0) = u(0);
    //     // g(0) = pinocchio::rnea(model, data, x.head(7), x.segment(7, 7), u);

    // }

    // template<typename T> 
    // EIGEN_STRONG_INLINE void
    // inequality_constraints_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
    //                             const Ref<const parameter_t <T>> p, const Ref<const static_parameter_t> d,
    //                             const scalar_t &t, Ref<constraint_t < T>>g) const noexcept
    // {
    //     g(0) = (T)0.0;
    // }


    // template<typename T> 
    // template <typename T,
    //           std::enable_if_t<std::is_same<T,scalar_t>::value, bool> = true>
    // EIGEN_STRONG_INLINE void
    // inequality_constraints_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
    //                             const Ref<const parameter_t <T>> p, const Ref<const static_parameter_t> d,
    //                             const scalar_t &t, Ref<constraint_t < T>>g) const noexcept
    // {
    //     // Matrix<double, 7, 1> q;
    //     // Matrix<double, 7, 1> q_dot;
    //     // Matrix<double, 7, 1> q_dot_dot;
    //     Eigen::Ref<Eigen::Matrix<double, 7, 1>> q(x.head(7));
    //     Eigen::Ref<Eigen::Matrix<double, 7, 1>> q_dot(x.tail(7));
    //     Eigen::Ref<Eigen::Matrix<double, 7, 1>> q_dot_dot(u);

    //     g(0) = pinocchio::rnea(model, data, q, q_dot, q_dot_dot);
    // }

    // // template<typename T> 
    // template <typename T,
    //           std::enable_if_t<std::is_same<T,ad_scalar_t>::value, bool> = true>
    // EIGEN_STRONG_INLINE void
    // inequality_constraints_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
    //                             const Ref<const parameter_t <T>> p, const Ref<const static_parameter_t> d,
    //                             const scalar_t &t, Ref<constraint_t < T>>g) const noexcept
    // {
    //     g(0) = (T)0.0;
    // }

    // // template<typename T> 
    // template <typename T,
    //           std::enable_if_t<std::is_same<T,ad2_scalar_t>::value, bool> = true>
    // EIGEN_STRONG_INLINE void
    // inequality_constraints_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
    //                             const Ref<const parameter_t <T>> p, const Ref<const static_parameter_t> d,
    //                             const scalar_t &t, Ref<constraint_t < T>>g) const noexcept
    // {
    //     g(0) = (T)0.0;
    // }

    EIGEN_STRONG_INLINE void
    inequality_constraints_impl(const Ref<const state_t<scalar_t>> x, const Ref<const control_t<scalar_t>> u,
                                const Ref<const parameter_t <scalar_t>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, Ref<constraint_t <scalar_t>>g) const noexcept
    {

        Eigen::Matrix<double, 7, 1> q = x.head(7);
        Eigen::Matrix<double, 7, 1> q_dot = x.tail(7);
        Eigen::Matrix<double, 7, 1> q_dot_dot = u;

        pinocchio::rnea(model, data, q, q_dot, q_dot_dot);

    }

    // template<typename T> 
    EIGEN_STRONG_INLINE void
    inequality_constraints_impl(const Ref<const state_t<ad_scalar_t>> x, const Ref<const control_t<ad_scalar_t>> u,
                                const Ref<const parameter_t <ad_scalar_t>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, Ref<constraint_t <ad_scalar_t>>g) const noexcept
    {
        // g(0) = (T)0.0;
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