#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    ConZono::ConZono(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                     const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                     const bool zero_one_form)
    {
        set(G, c, A, b, zero_one_form);
        sharp = true;
    }

    HybZono* ConZono::clone() const
    {
        return new ConZono(*this);
    }

    void ConZono::set(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                      const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                      const bool zero_one_form)
    {
        // check dimensions
        if (G.rows() != c.size() || A.rows() != b.size() || G.cols() != A.cols())
        {
            throw std::invalid_argument("ConZono: inconsistent dimensions.");
        }

        // conzono parameters
        this->G = G;
        this->A = A;
        this->c = c;
        this->b = b;
        this->nG = static_cast<int>(G.cols());
        this->nC = static_cast<int>(A.rows());
        this->n = static_cast<int>(G.rows());
        this->zero_one_form = zero_one_form;

        // abstract zono parameters
        this->nGc = this->nG;
        this->nGb = 0;
        this->Gc = this->G;
        this->Gb.resize(this->n, 0);
        this->Ac = this->A;
        this->Ab.resize(0, 0);
    }

    void ConZono::convert_form()
    {
        Eigen::Vector<zono_float, -1> c, b;
        Eigen::SparseMatrix<zono_float> G, A;

        if (!this->zero_one_form) // convert to [0,1] generators
        {
            c = this->c - this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            b = this->b + this->A * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            G = 2.0 * this->G;
            A = 2.0 * this->A;

            set(G, c, A, b, true);
        }
        else // convert to [-1,1] generators
        {
            c = this->c + 0.5 * this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            b = this->b - 0.5 * this->A * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            G = 0.5 * this->G;
            A = 0.5 * this->A;

            set(G, c, A, b, false);
        }
    }

    std::string ConZono::print() const
    {
        std::stringstream ss;
        ss << "ConZono: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "nG: " << this->nG << std::endl;
        ss << "nC: " << this->nC << std::endl;
        ss << "G: " << Eigen::Matrix<zono_float, -1, -1>(this->G) << std::endl;
        ss << "c: " << this->c << std::endl;
        ss << "A: " << Eigen::Matrix<zono_float, -1, -1>(this->A) << std::endl;
        ss << "b: " << this->b << std::endl;
        ss << "zero_one_form: " << this->zero_one_form;
        return ss.str();
    }

    Eigen::Vector<zono_float, -1> ConZono::do_optimize_over(
        const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q, const zono_float c,
        const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
        const WarmStartParams& warm_start_params) const
    {
        // check dimensions
        if (P.rows() != this->n || P.cols() != this->n || q.size() != this->n)
        {
            throw std::invalid_argument("Optimize over: inconsistent dimensions.");
        }

        // cost
        Eigen::SparseMatrix<zono_float> P_fact = this->G.transpose() * P * this->G;
        Eigen::Vector<zono_float, -1> q_fact = this->G.transpose() * (P * this->c + q);
        zono_float delta_c = (0.5 * this->c.transpose() * P * this->c + q.transpose() * this->c)(0);

        // solve QP
        OptSolution sol = this->qp_opt(P_fact, q_fact, c + delta_c, this->A, this->b, settings, solution,
                                       warm_start_params);

        // check feasibility and return solution
        if (sol.infeasible)
            return Eigen::Vector<zono_float, -1>(this->nG);
        else
            return this->G * sol.z + this->c;
    }

    Eigen::Vector<zono_float, -1> ConZono::do_project_point(const Eigen::Vector<zono_float, -1>& x,
                                                            const OptSettings& settings,
                                                            std::shared_ptr<OptSolution>* solution,
                                                            const WarmStartParams& warm_start_params) const
    {
        // check dimensions
        if (this->n != x.size())
        {
            throw std::invalid_argument("Point projection: inconsistent dimensions.");
        }

        // cost
        Eigen::SparseMatrix<zono_float> P = this->G.transpose() * this->G;
        Eigen::Vector<zono_float, -1> q = this->G.transpose() * (this->c - x);

        // solve QP
        const OptSolution sol = this->qp_opt(P, q, 0, this->A, this->b, settings, solution, warm_start_params);

        // check feasibility and return solution
        if (sol.infeasible)
            throw std::invalid_argument("Point projection: infeasible");
        else
            return this->G * sol.z + this->c;
    }

    bool ConZono::do_is_empty(const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                              const WarmStartParams& warm_start_params) const
    {
        // trivial case
        if (this->n == 0)
            return true;

        // cost
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        P.setIdentity();
        Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(this->nG);

        // solve QP
        OptSolution sol = this->qp_opt(P, q, 0, this->A, this->b, settings, solution, warm_start_params);

        // check infeasibility flag
        return sol.infeasible;
    }

    zono_float ConZono::do_support(const Eigen::Vector<zono_float, -1>& d,
                                   const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                   const WarmStartParams& warm_start_params)
    {
        // check dimensions
        if (this->n != d.size())
        {
            throw std::invalid_argument("Support: inconsistent dimensions.");
        }

        // cost
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // solve QP
        const OptSolution sol = this->qp_opt(P, q, 0, this->A, this->b, settings, solution, warm_start_params);

        // check feasibility and return solution
        if (sol.infeasible) // Z is empty
            throw std::invalid_argument("Support: infeasible");
        else
            return d.dot(this->G * sol.z + this->c);
    }

    bool ConZono::do_contains_point(const Eigen::Vector<zono_float, -1>& x,
                                    const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                    const WarmStartParams& warm_start_params) const
    {
        // check dimensions
        if (this->n != x.size())
        {
            throw std::invalid_argument("Contains point: inconsistent dimensions.");
        }

        // build QP for ADMM
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        P.setIdentity();
        Eigen::Vector<zono_float, -1> q(this->nG);
        q.setZero(); // zeros
        Eigen::SparseMatrix<zono_float> A = vcat<zono_float>(this->A, this->G);
        Eigen::Vector<zono_float, -1> b(this->nC + this->n);
        b.segment(0, this->nC) = this->b;
        b.segment(this->nC, this->n) = x - this->c;

        const OptSolution sol = this->qp_opt(P, q, 0, A, b, settings, solution, warm_start_params);

        // check feasibility and return solution
        return !(sol.infeasible);
    }

    OptSolution ConZono::qp_opt(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                                const zono_float c, const Eigen::SparseMatrix<zono_float>& A,
                                const Eigen::Vector<zono_float, -1>& b,
                                const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                const WarmStartParams& warm_start_params) const
    {
        // setup QP
        Eigen::Vector<zono_float, -1> xi_lb(this->nG);
        if (this->zero_one_form)
            xi_lb.setZero();
        else
            xi_lb.setConstant(-1);
        const Eigen::Vector<zono_float, -1> xi_ub = Eigen::Vector<zono_float, -1>::Ones(this->nG);

        const auto data = std::make_shared<ADMM_data>(P, q, A, b, xi_lb, xi_ub, c, settings);
        ADMM_solver solver(data);

        // warm start if applicable
        if (warm_start_params.z.size() == this->nG)
        {
            const Eigen::Vector<zono_float, -1> u0 = warm_start_params.u.size() == this->nG
                                                         ? warm_start_params.u
                                                         : Eigen::Vector<zono_float, -1>::Zero(this->nG);
            solver.warmstart(warm_start_params.z, u0);
        }

        // solve
        OptSolution sol = solver.solve();
        if (solution != nullptr)
            *solution = std::make_shared<OptSolution>(sol);
        return sol;
    }

    // bounding box
    Box ConZono::do_bounding_box(const OptSettings& settings, std::shared_ptr<OptSolution>*,
                                 const WarmStartParams& warm_start_params)
    {
        // make sure dimension is at least 1
        if (this->n == 0)
        {
            throw std::invalid_argument("Bounding box: empty set");
        }

        // init search direction for bounding box
        Eigen::Vector<zono_float, -1> d(this->n);
        d.setZero();

        // declarations
        Box box(this->n); // init
        // declare
        zono_float s_neg, s_pos;

        // build QP for ADMM
        const Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // convex solution
        Eigen::Vector<zono_float, -1> xi_lb, xi_ub;
        if (this->zero_one_form)
            xi_lb = Eigen::Vector<zono_float, -1>::Zero(this->nG);
        else
            xi_lb = -1.0 * Eigen::Vector<zono_float, -1>::Ones(this->nG);

        xi_ub = Eigen::Vector<zono_float, -1>::Ones(this->nG);

        // build ADMM object
        const auto data = std::make_shared<ADMM_data>(P, q, this->A, this->b, xi_lb, xi_ub, zero, settings);
        ADMM_solver solver(data);

        // check if warm start is valid
        const bool warm_start_valid = warm_start_params.z.size() == this->nG && warm_start_params.u.size() == this->nG;

        // get support in all box directions
        for (int i = 0; i < this->n; i++)
        {
            // negative direction

            // update QP
            d.setZero();
            d(i) = -1;
            data->q = -this->G.transpose() * d;

            // solve
            if (warm_start_valid)
            {
                solver.warmstart(warm_start_params.z, warm_start_params.u);
            }
            OptSolution sol = solver.solve();
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_neg = -d.dot(this->G * sol.z + this->c);

            // positive direction

            // update QP
            d.setZero();
            d(i) = 1;
            data->q = -this->G.transpose() * d;

            // solve
            if (warm_start_valid)
            {
                solver.warmstart(warm_start_params.z, warm_start_params.u);
            }
            sol = solver.solve();
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_pos = d.dot(this->G * sol.z + this->c);

            // store bounds
            box[i] = Interval(s_neg, s_pos);
        }

        return box;
    }
}
