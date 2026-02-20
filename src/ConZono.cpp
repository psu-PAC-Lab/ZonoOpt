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

    std::unique_ptr<Zono> ConZono::to_zono_approx() const
    {
        // check for case that there are no constraints
        if (this->nG == 0)
        {
            return std::make_unique<Point>(this->c);
        }
        if (this->nC == 0)
        {
            return std::make_unique<Zono>(this->G, this->c, this->zero_one_form);
        }

        // compute SVD of A
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        const Eigen::BDCSVD<Eigen::Matrix<zono_float, -1, -1>, Eigen::ComputeFullV | Eigen::ComputeFullU> svd(
            this->A.toDense());
#else
        const Eigen::BDCSVD<Eigen::Matrix<zono_float, -1, -1>> svd(
            this->A.toDense(), Eigen::ComputeFullV | Eigen::ComputeFullU);
#endif
        const Eigen::Vector<zono_float, -1>& sin_vals = svd.singularValues();

        const int n = this->nG;
        const int r = static_cast<int>(svd.rank());
        const int d = n - r;
        const Eigen::Matrix<zono_float, -1, -1> Vr = svd.matrixV().block(0, 0, n, r);
        const Eigen::Matrix<zono_float, -1, -1> Vd = svd.matrixV().block(0, r, n, d);

        // get xi_tilde_r
        const Eigen::Vector<zono_float, -1> b_tilde = svd.matrixU().transpose() * this->b;
        const Eigen::Vector<zono_float, -1> xi_tilde_r = b_tilde.segment(0, r).array() / sin_vals.segment(0, r).array();

        // build bounding zonotope
        const Eigen::SparseMatrix<zono_float> Vd_VdT_sp = (Vd * Vd.transpose()).sparseView();
        const Eigen::SparseMatrix<zono_float> G_zono = this->G * Vd_VdT_sp;
        const Eigen::Vector<zono_float, -1> c_zono = this->c + this->G * (Vr * xi_tilde_r);
        return std::make_unique<Zono>(G_zono, c_zono, this->zero_one_form);
    }

    void ConZono::constraint_reduction()
    {
        // make sure there are constraints to remove
        if (this->nC == 0) return;

        // put set into [-1, 1] form
        if (this->zero_one_form) this->convert_form();

        // execute algorithm 1 from paper
        Eigen::Vector<zono_float, -1> x_lb(this->nG);
        Eigen::Vector<zono_float, -1> x_ub(this->nG);
        x_lb.setConstant(-1);
        x_ub.setConstant(1);
        Box E(x_lb, x_ub);
        x_lb.setConstant(-std::numeric_limits<zono_float>::infinity());
        x_ub.setConstant(std::numeric_limits<zono_float>::infinity());
        Box R(x_lb, x_ub);
        Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_rm = this->A;
        for (int i = 0; i < this->nC; ++i)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_j(A_rm, i); it_j; ++it_j)
            {
                const zono_float a_ij = it_j.value();
                Interval y(this->b(i) / a_ij, this->b(i) / a_ij);
                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_k(A_rm, i); it_k; ++it_k)
                {
                    if (it_j.col() == it_k.col()) continue;
                    y = y - E[it_k.col()].to_interval() * (it_k.value() / a_ij);
                }
                R[it_j.col()].intersect_assign(R[it_j.col()], y.as_view());
                E[it_j.col()].intersect_assign(E[it_j.col()], R[it_j.col()]);
            }
        }

        // make sure conzono isn't empty (interval check)
        for (int j = 0; j < this->nG; ++j)
        {
            if (E[j].is_empty())
                throw std::runtime_error("ConZono constraint reduction: set is empty");
        }

        // build Q matrix
        Eigen::SparseMatrix<zono_float> Q(this->nG + this->nC, this->nG + this->nC);

        Eigen::SparseMatrix<zono_float> I_nG(this->nG, this->nG);
        I_nG.setIdentity();
        const Eigen::SparseMatrix<zono_float> Phi = this->G.transpose() * this->G + I_nG;

        std::vector<Eigen::Triplet<zono_float>> triplets;
        get_triplets_offset<zono_float>(Phi, triplets, 0, 0);
        get_triplets_offset<zono_float>(this->A, triplets, this->nG, 0);
        get_triplets_offset<zono_float>(this->A.transpose(), triplets, 0, this->nG);
        Q.setFromTriplets(triplets.begin(), triplets.end());

        // factorize Q
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>> Q_ldlt(Q);
        if (Q_ldlt.info() != Eigen::Success)
            throw std::runtime_error(
                "ConZono constraint reduction: Q matrix factorization failed, most likely A is not full row rank.");

        // get estimated Hausdorff error for eliminating each generator
        std::vector<std::pair<int, zono_float>> haus_vec; // (index, error)
        haus_vec.reserve(this->nG);
        auto shift_permute = [size=this->nG + this->nC + 1](const int start_index,
                                                            const int end_index) -> Eigen::PermutationMatrix<
            Eigen::Dynamic, Eigen::Dynamic>
        {
            assert(
                start_index >= 0 && end_index >= 0 && start_index < size && end_index < size && start_index <
                end_index);
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(size);
            for (int i = 0; i < start_index; ++i)
            {
                P.indices()[i] = i;
            }
            P.indices()[start_index] = end_index;
            for (int i = start_index + 1; i <= end_index; ++i)
            {
                P.indices()[i] = i - 1;
            }
            for (int i = end_index + 1; i < size; ++i)
            {
                P.indices()[i] = i;
            }
            return P;
        };
        Eigen::Vector<zono_float, -1> e_j(this->nG + this->nC);
        for (int j = 0; j < this->nG; ++j)
        {
            // get r_j
            const zono_float r_j = std::max<zono_float>(
                zero, std::max<zono_float>(std::abs(R[j].to_interval().lb), std::abs(R[j].to_interval().ub)) - one);
            if (r_j < zono_eps)
            {
                haus_vec.emplace_back(j, zero);
                continue;
            }

            // linear system matrix
            e_j.setZero();
            e_j(j) = 1;
            const Eigen::Vector<zono_float, -1> Qinv_e_j = Q_ldlt.solve(e_j);

            triplets.clear();
            for (int i = 0; i < this->nG + this->nC; ++i)
            {
                triplets.emplace_back(i, i, one);
                if (i == j)
                {
                    triplets.emplace_back(this->nG + this->nC, j, one); // extra row for e_j^T
                }
            }
            for (int i = 0; i < this->nG + this->nC; ++i)
            {
                triplets.emplace_back(i, this->nG + this->nC, Qinv_e_j(i));
            }
            Eigen::SparseMatrix<zono_float> M(this->nG + this->nC + 1, this->nG + this->nC + 1);
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
            M.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
            M.setFromTriplets(triplets.begin(), triplets.end());
#endif

            // RHS for linear system
            Eigen::Vector<zono_float, -1> rhs(this->nG + this->nC + 1);
            rhs.setZero();
            rhs(this->nG + this->nC) = r_j;

            // permutation matrices to make system upper triangular
            const auto P_R = shift_permute(j, this->nG + this->nC);
            const auto P_L = shift_permute(j, this->nG + this->nC - 1);

            // solve linear system in permuted space
            const Eigen::SparseMatrix<zono_float> M_perm = P_L * M * P_R.inverse();
            const Eigen::Vector<zono_float, -1> rhs_perm = P_L * rhs;
            const Eigen::Vector<zono_float, -1> y_perm = M_perm.triangularView<Eigen::Upper>().solve(rhs_perm);
            const Eigen::Vector<zono_float, -1> y = P_R * y_perm;

            // get Hausdorff distance estimate
            const Eigen::Vector<zono_float, -1> d = y.segment(0, this->nG);
            const zono_float haus_j = (this->G * d).norm() + d.norm();
            haus_vec.emplace_back(j, haus_j);
        }

        // sort by Hausdorff error
        std::sort(haus_vec.begin(), haus_vec.end(), [](const auto& a, const auto& b) { return a.second < b.second; });

        Eigen::SparseMatrix<zono_float> Ea(this->nG, this->nC); // init
        int gen_remove = -1;
        int cons_remove = -1;
        for (const auto& [j, err] : haus_vec)
        {
            // loop through column k of A to find a constraint to remove
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->A, j); it; ++it)
            {
                if (std::abs(it.value()) > zono_eps)
                {
                    triplets.emplace_back(j, static_cast<int>(it.row()), one / it.value());
                    Ea.insert(j, it.row()) = 1 / it.value();
                    gen_remove = j;
                    cons_remove = static_cast<int>(it.row());
                    break;
                }
            }

            // check if done
            if (gen_remove != -1 && cons_remove != -1)
                break;
        }
        if (gen_remove == -1 || cons_remove == -1)
            throw std::runtime_error("ConZono: constraint reduction cannot find valid constraint to remove");

        // apply algorithm from Scott paper
        const Eigen::SparseMatrix<zono_float> Lambda_G = G * Ea;
        const Eigen::SparseMatrix<zono_float> Lambda_A = A * Ea;
        Eigen::SparseMatrix<zono_float> Gp = this->G - Lambda_G * this->A;
        Eigen::Vector<zono_float, -1> cp = this->c + Lambda_G * this->b;
        Eigen::SparseMatrix<zono_float> Ap = this->A - Lambda_A * this->A;
        Eigen::Vector<zono_float, -1> bp = this->b - Lambda_A * this->b;

        // generator removal matrix
        triplets.clear();
        for (int j = 0; j < gen_remove; ++j)
        {
            triplets.emplace_back(j, j, one);
        }
        for (int j = gen_remove + 1; j < this->nG; ++j)
        {
            triplets.emplace_back(j, j - 1, one);
        }
        Eigen::SparseMatrix<zono_float> dG(this->nG, this->nG - 1);
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        dG.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        dG.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // constraint removal matrix
        triplets.clear();
        for (int i = 0; i < cons_remove; ++i)
        {
            triplets.emplace_back(i, i, one);
        }
        for (int i = cons_remove + 1; i < this->nC; ++i)
        {
            triplets.emplace_back(i - 1, i, one);
        }
        Eigen::SparseMatrix<zono_float> dA(this->nC - 1, this->nC);
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        dA.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        dA.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // update
        this->set(Gp * dG, cp, dA * Ap * dG, dA * bp, false);
    }

    std::unique_ptr<ConZono> vrep_2_conzono(const Eigen::Matrix<zono_float, -1, -1>& Vpoly)
    {
        // dimensions
        const int n_dims = static_cast<int>(Vpoly.cols());
        const int n_verts = static_cast<int>(Vpoly.rows());

        // make generators
        const Eigen::SparseMatrix<zono_float> G = Vpoly.transpose().sparseView();
        Eigen::Vector<zono_float, -1> c(n_dims);
        c.setZero();

        // make constraints
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        Eigen::SparseMatrix<zono_float> A(1, n_verts);
        for (int i = 0; i < n_verts; i++)
        {
            tripvec.emplace_back(0, i, one);
        }
        A.setFromTriplets(tripvec.begin(), tripvec.end());

        Eigen::Vector<zono_float, -1> b(1);
        b(0) = 1;

        // return conzono
        return std::make_unique<ConZono>(G, c, A, b, true);
    }

    std::unique_ptr<HybZono> ConZono::do_complement(const zono_float delta_m, bool, const OptSettings&,
                                                    std::shared_ptr<OptSolution>*, int, int)
    {
        // make sure in [-1,1] form
        if (this->is_0_1_form()) this->convert_form();

        // get a value lambda_m such that lambda_m >= max{ ||lambda||_infty : |[G^T A^T] lambda| <= 1 }

        // construct the matrix [G^T A^T]
        const auto GT = this->G.transpose();
        const auto AT = this->A.transpose();
        const Eigen::SparseMatrix<zono_float, Eigen::RowMajor> GTAT = hcat<zono_float>(GT, AT); // convert to row-major

        // get the smallest non-zero element in the matrix, making sure that all rows have at least one non-zero element
        zono_float min_non_zero = std::numeric_limits<zono_float>::max();
        for (int row = 0; row < GTAT.outerSize(); ++row)
        {
            bool value_exists = false;
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(GTAT, row); it; ++it)
            {
                if (it.value() < min_non_zero && std::abs(it.value()) > zono_eps)
                {
                    min_non_zero = it.value();
                }
                value_exists = true;
            }
            if (!value_exists)
            {
                std::stringstream ss;
                ss << "ConZono complement: row " << row << " of [G^T A^T] has no non-zero elements.";
                throw std::runtime_error(ss.str());
            }
        }

        // value for lambda_m
        const zono_float lambda_m = 1 / min_non_zero;

        // m value
        const zono_float m = delta_m + 1;

        // interval sets
        Eigen::Vector<zono_float, -1> l1(2 * this->nG);
        l1.setConstant(-(m + delta_m / 2));
        Eigen::Vector<zono_float, -1> u1(2 * this->nG);
        u1.setConstant(1 + delta_m / 2);
        auto Zf1 = interval_2_zono(Box(l1, u1));
        if (Zf1->is_0_1_form()) Zf1->convert_form();

        Eigen::Vector<zono_float, -1> l2(4 * this->nG);
        l2.segment(0, 2 * this->nG).setConstant(-(m + 3 * delta_m / 2 + 1));
        l2.segment(2 * this->nG, 2 * this->nG).setConstant(-2);
        Eigen::Vector<zono_float, -1> u2(4 * this->nG);
        u2.segment(0, 2 * this->nG).setConstant(delta_m / 2);
        u2.segment(2 * this->nG, 2 * this->nG).setZero();
        auto Zf2 = interval_2_zono(Box(l2, u2));
        if (Zf2->is_0_1_form()) Zf2->convert_form();

        // build complement

        // Gc
        Eigen::SparseMatrix<zono_float> Gc = m * this->G;
        Gc.conservativeResize(this->n, 9 * this->nG + this->n + this->nC + 1);

        // Gb
        Eigen::SparseMatrix<zono_float> Gb(this->n, 2 * this->nG);

        // c
        Eigen::Vector<zono_float, -1> c = this->c;

        // helper matrices

        // declarations
        std::vector<Eigen::Triplet<zono_float>> triplets;
        int n_offset = 0;

        // AcPF = [mI, -delta_m/2, 0, 0, 0;
        //         -mI, -delta_m/2, 0, 0, 0]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, this->nG, -delta_m / two);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, i, -m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG, -delta_m / two);
        }
        Eigen::SparseMatrix<zono_float> AcPF(2 * this->nG, 3 * this->nG + 1 + this->nC + this->n);
        AcPF.setFromTriplets(triplets.begin(), triplets.end());

        // AcDF = [0, 0, lambda_m [G^T A^T], 0.5 I, -0.5 I;
        //         0, 0, 0, 0.5 1^T, 0.5 1^T]
        triplets.clear();
        get_triplets_offset<zono_float>(lambda_m * GTAT, triplets, 0, this->nG + 1);
        n_offset = this->nG + 1 + this->n + this->nC;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, n_offset + i, p5);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG, n_offset + i, p5);
        }
        n_offset += this->nG;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, n_offset + i, -p5);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG, n_offset + i, p5);
        }
        Eigen::SparseMatrix<zono_float> AcDF(this->nG + 1, 3 * this->nG + 1 + this->nC + this->n);
        AcDF.setFromTriplets(triplets.begin(), triplets.end());

        // AcCS = [-mI, delta_m/2, 0, 0, 0;
        //         mI, delta_m/2, 0, 0, 0;
        //         0, 0, 0, I, 0;
        //         0, 0, 0, 0, I]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, -m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, this->nG, delta_m / two);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG, delta_m / two);
        }
        n_offset = this->nG + 1 + this->n + this->nC;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(2 * this->nG + i, n_offset + i, one);
        }
        n_offset += this->nG;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(3 * this->nG + i, n_offset + i, one);
        }
        Eigen::SparseMatrix<zono_float> AcCS(4 * this->nG, 3 * this->nG + 1 + this->nC + this->n);
        AcCS.setFromTriplets(triplets.begin(), triplets.end());

        // AbCS = [mI, 0;
        //         0, mI;
        //         -I, 0;
        //         0, -I]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG + i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(2 * this->nG + i, i, -one);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(3 * this->nG + i, this->nG + i, -one);
        }
        Eigen::SparseMatrix<zono_float> AbCS(4 * this->nG, 2 * this->nG);
        AbCS.setFromTriplets(triplets.begin(), triplets.end());

        // bDF = [0;
        //        1-nG]
        Eigen::Vector<zono_float, -1> bDF(this->nG + 1);
        bDF.setZero();
        bDF(this->nG) = static_cast<zono_float>(1 - this->nG);

        // Ac = [mA, 0, 0, 0;
        //       AcPF, Gf1, 0;
        //       AcDF, 0, 0;
        //       AcCS, 0, Gf2]
        triplets.clear();
        int m_offset = 0;
        n_offset = this->nG + 1 + this->n + this->nC + 2 * this->nG;
        get_triplets_offset<zono_float>(m * this->A, triplets, 0, 0);
        m_offset += this->nC;
        get_triplets_offset<zono_float>(AcPF, triplets, m_offset, 0);
        get_triplets_offset<zono_float>(Zf1->G, triplets, m_offset, n_offset);
        m_offset += 2 * this->nG;
        get_triplets_offset<zono_float>(AcDF, triplets, m_offset, 0);
        m_offset += this->nG + 1;
        get_triplets_offset<zono_float>(AcCS, triplets, m_offset, 0);
        get_triplets_offset<zono_float>(Zf2->G, triplets, m_offset, n_offset + 2 * this->nG);
        Eigen::SparseMatrix<zono_float> Ac(7 * this->nG + this->nC + 1, 9 * this->nG + this->n + this->nC + 1);
        Ac.setFromTriplets(triplets.begin(), triplets.end());

        // Ab = [0;
        //       0;
        //       0;
        //       AbCS];
        triplets.clear();
        get_triplets_offset<zono_float>(AbCS, triplets, this->nC + 2 * this->nG + this->nG + 1, 0);
        Eigen::SparseMatrix<zono_float> Ab(7 * this->nG + this->nC + 1, 2 * this->nG);
        Ab.setFromTriplets(triplets.begin(), triplets.end());

        // b = [b;
        //      cf1;
        //      bDF;
        //      cf2]
        Eigen::Vector<zono_float, -1> b(7 * this->nG + this->nC + 1);
        b.segment(0, this->nC) = this->b;
        b.segment(this->nC, 2 * this->nG) = Zf1->c;
        b.segment(this->nC + 2 * this->nG, this->nG + 1) = bDF;
        b.segment(this->nC + 2 * this->nG + this->nG + 1, 4 * this->nG) = Zf2->c;

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, false, false);
    }


}
