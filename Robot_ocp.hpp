#include "polynomials/ebyshev.hpp"
#include "polynomials/splines.hpp"
#include "control/continuous_ocp.hpp"
#include "utils/helpers.hpp"

#include "solvers/sqp_base.hpp"
#include "solvers/osqp_interface.hpp"
#include "control/mpc_wrapper.hpp"

#include <iomanip>
#include <iostream>
#include <chrono>

#define POLY_ORDER 5
#define NUM_SEG    2

/** benchmark the new collocation class */
using Polynomial = polympc::Chebyshev<POLY_ORDER, polympc::GAUSS_LOBATTO, double>;
using Approximation = polympc::Spline<Polynomial, NUM_SEG>;

POLYMPC_FORWARD_DECLARATION(/*Name*/ RobotOCP, /*NX*/ 3, /*NU*/ 2, /*NP*/ 0, /*ND*/ 1, /*NG*/0, /*TYPE*/ double)

using namespace Eigen;

class RobotOCP : public ContinuousOCP<RobotOCP, Approximation, SPARSE>
{
public:
    ~RobotOCP() = default;

    RobotOCP()
    {
        /** initialise weight matrices to identity (for example)*/
        Q.setIdentity();
        R.setIdentity();
        QN.setIdentity();

    }

    Matrix<scalar_t, 3, 3> Q;
    Matrix<scalar_t, 2, 2> R;
    Matrix<scalar_t, 3, 3> QN;

    template<typename T>
    inline void dynamics_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
                              const Ref<const parameter_t<T>> p, const Ref<const static_parameter_t> &d,
                              const T &t, Ref<state_t<T>> xdot) const noexcept
    {
        polympc::ignore_unused_var(p);
        polympc::ignore_unused_var(t);

        xdot(0) = u(0) * cos(x(2)) * cos(u(1));
        xdot(1) = u(0) * sin(x(2)) * cos(u(1));
        xdot(2) = u(0) * sin(u(1)) / d(0);
    }

    template<typename T>
    inline void lagrange_term_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
                                   const Ref<const parameter_t<T>> p, const Ref<const static_parameter_t> d,
                                   const scalar_t &t, T &lagrange) noexcept
    {
        polympc::ignore_unused_var(p);
        polympc::ignore_unused_var(t);
        polympc::ignore_unused_var(d);

        lagrange = x.dot(Q.template cast<T>() * x) + u.dot(R.template cast<T>() * u);
    }

    template<typename T>
    inline void mayer_term_impl(const Ref<const state_t<T>> x, const Ref<const control_t<T>> u,
                                const Ref<const parameter_t<T>> p, const Ref<const static_parameter_t> d,
                                const scalar_t &t, T &mayer) noexcept
    {
        polympc::ignore_unused_var(p);
        polympc::ignore_unused_var(t);
        polympc::ignore_unused_var(d);
        polympc::ignore_unused_var(u);

        mayer = x.dot(QN.template cast<T>() * x);
    }
};