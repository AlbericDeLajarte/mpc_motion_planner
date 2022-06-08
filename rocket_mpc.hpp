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

#define POLY_ORDER 5
#define NUM_SEG    2
#define NUM_EXP    1

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ guidance_ocp, /*NX*/ 8, /*NU*/ 1, /*NP*/ 0, /*ND*/ 0, /*NG*/0, /*TYPE*/ double)
using namespace Eigen;

class guidance_ocp : public ContinuousOCP<guidance_ocp, Approximation, DENSE>
{
public:
    ~guidance_ocp() = default;

    double cos_vertical_angle;

    template<typename T>
    inline void dynamics_impl(const Eigen::Ref<const state_t<T>> x, const Eigen::Ref<const control_t<T>> u,
                              const Eigen::Ref<const parameter_t<T>> p, const Eigen::Ref<const static_parameter_t> &d,
                              const T &t, Eigen::Ref<state_t<T>> xdot) const noexcept
    {
        // -------------- Simulation variables -----------------------------
        T mass = 20 + x(6);                   // Instantaneous mass of the rocekt in [kg]
        Matrix<T, 3, 1> gravity;
        gravity <<   (T)0.0, (T)0.0, (T)9.81;

        Matrix<T, 3, 1> thrust_force; thrust_force << (T)0.0, (T)0.0, u(0);

        // -------------- Differential equation ---------------------

        // Position variation is mass
        xdot.head(3) = x.segment(3,3);

        // Speed variation is Force/mass
        xdot.segment(3,3) = thrust_force/mass - gravity;

        // Mass variation is proportional to thrust
        xdot(6) = -u(0)/(200*gravity[2]);  // When parameters are working -> Use if(t<p(0))

        xdot(7) = 0;
                
        xdot *= x(7);
        
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
        mayer = x(7);

        // polympc::ignore_unused_var(x);
        polympc::ignore_unused_var(u);
        polympc::ignore_unused_var(d);
        polympc::ignore_unused_var(t);
  
    }
};

#endif //SRC_ROCKET_MPC_HPP