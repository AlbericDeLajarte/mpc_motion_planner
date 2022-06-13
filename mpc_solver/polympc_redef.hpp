#ifndef SRC_POLYMPC_REDEF_HPP
#define SRC_POLYMPC_REDEF_HPP

#include "polynomials/splines.hpp"

#include "solvers/sqp_base.hpp"
#include "solvers/box_admm.hpp"
//#include "solvers/osqp_interface.hpp"

#include "utils/helpers.hpp"

using namespace Eigen;


// ------ Create solver -------
template<typename Problem, typename QPSolver>
class MySolver;

template<typename Problem, typename QPSolver = boxADMM<Problem::VAR_SIZE, Problem::DUAL_SIZE, typename Problem::scalar_t>>
class MySolver : public SQPBase<MySolver<Problem, QPSolver>, Problem, QPSolver> {
public:
    using Base = SQPBase<MySolver<Problem, QPSolver>, Problem, QPSolver>;
    using typename Base::scalar_t;
    using typename Base::nlp_variable_t;
    using typename Base::nlp_hessian_t;

    // using typename Base::scalar_t;
    using typename Base::parameter_t;
    // using typename Base::nlp_variable_t;
    using typename Base::nlp_dual_t;
    using typename Base::nlp_jacobian_t;
    using typename Base::nlp_constraints_t;
    // using typename Base::nlp_hessian_t;

    EIGEN_STRONG_INLINE Problem
    &

    get_problem() noexcept { return this->problem; }

    /** change Hessian regularization to non-default*/
    EIGEN_STRONG_INLINE void hessian_regularisation_dense_impl(
            Eigen::Ref<nlp_hessian_t> lag_hessian) noexcept { //Three methods: adding a constant identity matrix, estimating eigen values with Greshgoring circles, Eigen Value Decomposition
        const int n = this->m_H.rows();

        /**Regularize by the estimation of the minimum negative eigen value--does not work with inexact Hessian update(matrix is already PSD)*/
        scalar_t aii, ri;
        for (int i = 0; i < n; i++) {
            aii = lag_hessian(i, i);
            ri = (lag_hessian.col(i).cwiseAbs()).sum() -
                 abs(aii); // The hessian is symmetric, Greshgorin discs from rows or columns are equal
            if (aii - ri <= 0) {
                lag_hessian(i, i) += (ri - aii) + scalar_t(0.01);
            } //All Greshgorin discs are in the positive half
        }
    }

    EIGEN_STRONG_INLINE void hessian_regularisation_sparse_impl(
            nlp_hessian_t &lag_hessian) noexcept {//Three methods: adding a constant identity matrix, estimating eigen values with Gershgorin circles, Eigen Value Decomposition
        const int n = this->m_H.rows(); //132=m_H.toDense().rows()

        /**Regularize by the estimation of the minimum negative eigen value*/
        scalar_t aii, ri;
        for (int i = 0; i < n; i++) {
            aii = lag_hessian.coeffRef(i, i);
            ri = (lag_hessian.col(i).cwiseAbs()).sum() -
                 abs(aii); // The hessian is symmetric, Greshgorin discs from rows or columns are equal
            if (aii - ri <= 0)
                lag_hessian.coeffRef(i, i) += (ri - aii) + 0.001;//All Gershgorin discs are in the positive half
        }
    }

    /** change step size selection algorithm */
    scalar_t step_size_selection_impl(const Ref<const nlp_variable_t> &p) noexcept {
        //std::cout << "taking NEW implementation \n";
        scalar_t mu, phi_l1, Dp_phi_l1;
        nlp_variable_t cost_gradient = this->m_h;
        const scalar_t tau = this->m_settings.tau; // line search step decrease, 0 < tau < settings.tau

        scalar_t constr_l1 = this->constraints_violation(this->m_x);

        // TODO: get mu from merit function model using hessian of Lagrangian
        //const scalar_t quad_term = p.dot(this->m_H * p);
        //const scalar_t qt = quad_term >= 0 ? scalar_t(0.5) * quad_term : 0;
        //mu = (abs(cost_gradient.dot(p)) ) / ((1 - this->m_settings.rho) * constr_l1);

        mu = this->m_lam_k.template lpNorm<Eigen::Infinity>();

        //std::cout << "mu: " << mu << "\n";

        scalar_t cost_1;
        this->problem.cost(this->m_x, this->m_p, cost_1);

        //std::cout << "l1: " << constr_l1 << " cost: " << cost_1 << "\n";

        phi_l1 = cost_1 + mu * constr_l1;
        Dp_phi_l1 = cost_gradient.dot(p) - mu * constr_l1;

        scalar_t alpha = scalar_t(1.0);
        scalar_t cost_step;
        nlp_variable_t x_step;
        for (int i = 1; i < this->m_settings.line_search_max_iter; i++) {
            x_step.noalias() = alpha * p;
            x_step += this->m_x;
            this->problem.cost(x_step, this->m_p, cost_step);

            //std::cout << "i: " << i << " l1: " << this->constraints_violation(x_step) << " cost: " << cost_step << "\n";

            scalar_t phi_l1_step = cost_step + mu * this->constraints_violation(x_step);

            //std::cout << "phi before: " << phi_l1 << " after: " << phi_l1_step <<  " required diff: " << alpha * this->m_settings.eta * Dp_phi_l1 << "\n";

            if (phi_l1_step <= (phi_l1 + alpha * this->m_settings.eta * Dp_phi_l1)) {
                // accept step
                return alpha;
            } else {
                alpha = tau * alpha;
            }
        }

        return alpha;
    }

    /** change Hessian update algorithm to the one provided by ContinuousOCP*/
    EIGEN_STRONG_INLINE void
    hessian_update_impl(Eigen::Ref<nlp_hessian_t> hessian, const Eigen::Ref<const nlp_variable_t> &x_step,
                        const Eigen::Ref<const nlp_variable_t> &grad_step) noexcept {
        this->problem.hessian_update_impl(hessian, x_step, grad_step);
    }

    // ----------------- Fix from Roland ----------------- //
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv //

    EIGEN_STRONG_INLINE void
    update_linearisation_dense_impl(const Eigen::Ref<const nlp_variable_t>& x, const Eigen::Ref<const parameter_t>& p,
                                    const Eigen::Ref<const nlp_variable_t>& x_step, const Eigen::Ref<const nlp_dual_t>& lam,
                                    Eigen::Ref<nlp_variable_t> cost_grad, Eigen::Ref<nlp_hessian_t> lag_hessian,
                                    Eigen::Ref<nlp_jacobian_t> A, Eigen::Ref<nlp_constraints_t> b) noexcept{
        this->linearisation_dense_impl(x, p, lam, cost_grad, lag_hessian, A, b);
    }

    EIGEN_STRONG_INLINE void
    update_linearisation_sparse_impl(const Eigen::Ref<const nlp_variable_t>& x, const Eigen::Ref<const parameter_t>& p,
                                    const Eigen::Ref<const nlp_variable_t>& x_step, const Eigen::Ref<const nlp_dual_t>& lam,
                                    Eigen::Ref<nlp_variable_t> cost_grad, nlp_hessian_t& lag_hessian,
                                    nlp_jacobian_t& A, Eigen::Ref<nlp_constraints_t> b) noexcept{
        this->linearisation_sparse_impl(x, p, lam, cost_grad, lag_hessian, A, b);
    }

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ //
    // ----------------- Fix from Roland ----------------- //
    

};



#endif //SRC_POLYMPC_REDEF_HPP