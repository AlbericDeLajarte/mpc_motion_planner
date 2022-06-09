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
#define NUM_EXP    1

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ guidance_ocp, /*NX*/ 3, /*NU*/ 1, /*NP*/ 1, /*ND*/ 0, /*NG*/0, /*TYPE*/ double)
using namespace Eigen;

class guidance_ocp : public ContinuousOCP<guidance_ocp, Approximation, DENSE>
{
public:
    ~guidance_ocp() = default;

    template<typename T>
    inline void dynamics_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                              const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> &d,
                              const T &t, Eigen::Ref<state_t<T>> xdot) const noexcept
    {
        // -------------- Simulation variables -----------------------------
        T mass = (T)20;                   // Instantaneous mass of the rocekt in [kg]
        T g0 = (T)9.81;

        // -------------- Differential equation ---------------------

        // Position variation is mass
        xdot(0) = x(1);

        // Speed variation is Force/mass
        xdot(1) = u(0)/mass - g0;

        // Mass variation is proportional to thrust
        xdot(2) = -u(0)/(200*g0);  // When parameters are working -> Use if(t<p(0))
                
        xdot *= p(0);
        
        polympc::ignore_unused_var(t);
    }

    template<typename T>
    inline void lagrange_term_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                                   const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> d,
                                   const scalar_t &t, T &lagrange) noexcept
    {
        //lagrange = -1e-1*x(6);
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