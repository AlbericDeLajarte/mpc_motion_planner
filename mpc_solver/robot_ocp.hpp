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

// Global variable with rocket parameters and useful methods

using namespace std;

#define POLY_ORDER 10
#define NUM_SEG    1

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, float>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ guidance_ocp, /*NX*/ 14, /*NU*/ 7, /*NP*/ 0, /*ND*/ 0, /*NG*/0, /*TYPE*/ float)
using namespace Eigen;

class guidance_ocp : public ContinuousOCP<guidance_ocp, Approximation, DENSE>
{
public:
    ~guidance_ocp() = default;

    Matrix<scalar_t, 14, 14> Q;
    Matrix<scalar_t, 7, 7> R;
    Matrix<scalar_t, 14, 14> QN;

    Matrix<scalar_t, 14, 1> x_ref;

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
        // xdot *= p(0);
        
        polympc::ignore_unused_var(t);
    }

    template<typename T>
    inline void lagrange_term_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                                   const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> d,
                                   const scalar_t &t, T &lagrange) noexcept
    {
        Matrix<T, 14, 1> x_error = x - x_ref.template cast<T>();
        lagrange = x_error.dot(Q.template cast<T>() * x_error) + u.dot(R.template cast<T>() * u);
        
        // lagrange = (T)0.0;
    }

    template<typename T>
    inline void mayer_term_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                                const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> d,
                                const scalar_t &t, T &mayer) noexcept
    {   
        Matrix<T, 14, 1> x_error = x - x_ref.template cast<T>();
        mayer = x_error.dot(QN.template cast<T>() * x_error);

        // mayer = (T)0.0;

        // polympc::ignore_unused_var(x);
        // polympc::ignore_unused_var(u);
        // polympc::ignore_unused_var(d);
        // polympc::ignore_unused_var(t);
  
    }
};

#endif //SRC_ROCKET_MPC_HPP