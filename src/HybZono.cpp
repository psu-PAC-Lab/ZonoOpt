#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    HybZono::HybZono(const Eigen::SparseMatrix<zono_float>& Gc, const Eigen::SparseMatrix<zono_float>& Gb,
                     const Eigen::Vector<zono_float, -1>& c,
                     const Eigen::SparseMatrix<zono_float>& Ac, const Eigen::SparseMatrix<zono_float>& Ab,
                     const Eigen::Vector<zono_float, -1>& b,
                     const bool zero_one_form, const bool sharp)
    {
        set(Gc, Gb, c, Ac, Ab, b, zero_one_form, sharp);
    }

    HybZono* HybZono::clone() const
    {
        return new HybZono(*this);
    }

    void HybZono::set(const Eigen::SparseMatrix<zono_float>& Gc, const Eigen::SparseMatrix<zono_float>& Gb,
                      const Eigen::Vector<zono_float, -1>& c,
                      const Eigen::SparseMatrix<zono_float>& Ac, const Eigen::SparseMatrix<zono_float>& Ab,
                      const Eigen::Vector<zono_float, -1>& b,
                      const bool zero_one_form, const bool sharp)
    {
        // check dimensions
        if (Gc.rows() != c.size() || Gb.rows() != c.size() || Gc.cols() != Ac.cols()
            || Gb.cols() != Ab.cols() || Ac.rows() != b.size() || Ab.rows() != b.size())
        {
            throw std::invalid_argument("HybZono: inconsistent dimensions.");
        }

        this->Gc = Gc;
        this->Gb = Gb;
        this->Ac = Ac;
        this->Ab = Ab;
        this->c = c;
        this->b = b;
        this->nGc = static_cast<int>(Gc.cols());
        this->nGb = static_cast<int>(Gb.cols());
        this->nC = static_cast<int>(Ac.rows());
        this->n = static_cast<int>(Gc.rows());
        this->zero_one_form = zero_one_form;

        make_G_A();

        this->sharp = sharp;
    }

    void HybZono::convert_form()
    {
        Eigen::Vector<zono_float, -1> c, b;
        Eigen::SparseMatrix<zono_float> Gb, Ab, Ac, Gc;

        if (!this->zero_one_form) // convert to [0,1] generators
        {
            c = this->c - this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            b = this->b + this->A * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            Gb = 2.0 * this->Gb;
            Ab = 2.0 * this->Ab;
            Gc = 2.0 * this->Gc;
            Ac = 2.0 * this->Ac;

            set(Gc, Gb, c, Ac, Ab, b, true);
        }
        else // convert to [-1,1] generators
        {
            c = this->c + 0.5 * this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            b = this->b - 0.5 * this->A * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            Gb = 0.5 * this->Gb;
            Ab = 0.5 * this->Ab;
            Gc = 0.5 * this->Gc;
            Ac = 0.5 * this->Ac;

            set(Gc, Gb, c, Ac, Ab, b, false);
        }
    }

    bool HybZono::remove_redundancy(const int contractor_iter)
    {
        // initial complexity
        const int nG_init = this->nG;
        const int nC_init = this->nC;

        // declare vars
        std::set<int> idx_c_to_remove, idx_b_to_remove;

        // lambda to remove generators
        auto remove_all_generators = [this](const std::set<int>& idx_c, const std::set<int>& idx_b) -> void
        {
            // remove generators
            if (!idx_c.empty())
            {
                remove_generators(this->Gc, this->Ac, idx_c);
            }
            if (!idx_b.empty())
            {
                remove_generators(this->Gb, this->Ab, idx_b);
            }

            // update number of generators (needs to happen before call to make_G_A())
            this->nG = static_cast<int>(this->G.cols());
            this->nGc = static_cast<int>(this->Gc.cols());
            this->nGb = static_cast<int>(this->Gb.cols());

            // update equivalent matrices
            make_G_A();
        };

        // apply interval contractor
        Eigen::Vector<zono_float, -1> x_l(this->nG);
        Eigen::Vector<zono_float, -1> x_u(this->nG);
        if (this->zero_one_form)
        {
            x_l.setZero();
        }
        else
        {
            x_l.setConstant(-1);
        }
        x_u.setOnes();
        MI_Box box(x_l, x_u, {this->nGc, this->nGb}, this->zero_one_form);
        box.contract(this->A, this->b, contractor_iter);

        // find any variables whose values are fixed
        std::vector<std::pair<int, zono_float>> fixed_vars;
        for (int i = 0; i < this->nG; ++i)
        {
            if (box[i].is_single_valued())
            {
                fixed_vars.emplace_back(i, box[i].get_y_max());
                if (i < this->nGc)
                {
                    idx_c_to_remove.insert(i);
                }
                else
                {
                    idx_b_to_remove.insert(i - this->nGc);
                }
            }
        }

        // get updates to c and b
        Eigen::Vector<zono_float, -1> dc(this->n);
        Eigen::Vector<zono_float, -1> db(this->nC);
        Eigen::Vector<zono_float, -1> dc_k(this->n);
        Eigen::Vector<zono_float, -1> db_k(this->nC);
        dc.setZero();
        db.setZero();
        for (auto& [k, val] : fixed_vars)
        {
            dc_k.setZero();
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->G, k); it; ++it)
            {
                dc_k(it.row()) = it.value() * val;
            }
            dc += dc_k;

            db_k.setZero();
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->A, k); it; ++it)
            {
                db_k(it.row()) = it.value() * val;
            }
            db -= db_k;
        }

        // set center and constraint vector
        this->c += dc;
        this->b += db;

        // remove generators
        remove_all_generators(idx_c_to_remove, idx_b_to_remove);

        // remove redundant constraints
        remove_redundant_constraints<zono_float>(this->A, this->b);
        this->nC = static_cast<int>(this->A.rows());

        // update Ac, Ab
        set_Ac_Ab_from_A();

        // identify any unused generators
        idx_c_to_remove = find_unused_generators(this->Gc, this->Ac);
        idx_b_to_remove = find_unused_generators(this->Gb, this->Ab);

        // remove
        remove_all_generators(idx_c_to_remove, idx_b_to_remove);

        // flag indicating success
        return (this->nG < nG_init || this->nC < nC_init);
    }

    std::string HybZono::print() const
    {
        std::stringstream ss;
        ss << "HybZono: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "nGc: " << this->nGc << std::endl;
        ss << "nGb: " << this->nGb << std::endl;
        ss << "nC: " << this->nC << std::endl;
        ss << "Gc: " << Eigen::Matrix<zono_float, -1, -1>(this->Gc) << std::endl;
        ss << "Gb: " << Eigen::Matrix<zono_float, -1, -1>(this->Gb) << std::endl;
        ss << "c: " << this->c << std::endl;
        ss << "Ac: " << Eigen::Matrix<zono_float, -1, -1>(this->Ac) << std::endl;
        ss << "Ab: " << Eigen::Matrix<zono_float, -1, -1>(this->Ab) << std::endl;
        ss << "b: " << this->b << std::endl;
        ss << "zero_one_form: " << this->zero_one_form << std::endl;
        ss << "sharp: " << this->sharp;
        return ss.str();
    }

    std::ostream& operator<<(std::ostream& os, const HybZono& Z)
    {
        os << Z.print();
        return os;
    }

    Eigen::Vector<zono_float, -1> HybZono::do_optimize_over(
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

        // solve MIQP
        OptSolution sol = this->mi_opt(P_fact, q_fact, c + delta_c, this->A, this->b, settings, solution,
                                       warm_start_params);
        if (sol.infeasible)
            return Eigen::Vector<zono_float, -1>(this->nG);
        else
            return this->G * sol.z + this->c;
    }

    Eigen::Vector<zono_float, -1> HybZono::do_project_point(const Eigen::Vector<zono_float, -1>& x,
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

        // solve MIQP
        const OptSolution sol = this->mi_opt(P, q, 0, this->A, this->b, settings, solution, warm_start_params);
        if (sol.infeasible)
            throw std::runtime_error("Point projection: infeasible");

        return this->G * sol.z + this->c;
    }

    bool HybZono::do_is_empty(const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                              const WarmStartParams&) const
    {
        // trivial case
        if (this->n == 0)
            return true;

        // optimize over P=I, q=0
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        P.setIdentity();
        Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(this->nG);

        // solve
        const std::vector<OptSolution> sol_vec = this->
            mi_opt_multisol(P, q, 0, this->A, this->b, 1, settings, solution);
        if (sol_vec.size() > 0)
            return sol_vec[0].infeasible;
        else
            return true;
    }

    bool HybZono::do_contains_point(const Eigen::Vector<zono_float, -1>& x, const OptSettings& settings,
                                    std::shared_ptr<OptSolution>* solution,
                                    const WarmStartParams& warm_start_params) const
    {
        // check dimensions
        if (this->n != x.size())
        {
            throw std::invalid_argument("Contains point: inconsistent dimensions.");
        }

        // cost and constraints
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        P.setIdentity();
        Eigen::Vector<zono_float, -1> q(this->nG);
        q.setZero(); // zeros
        Eigen::SparseMatrix<zono_float> A = vcat<zono_float>(this->A, this->G);
        Eigen::Vector<zono_float, -1> b(this->nC + this->n);
        b.segment(0, this->nC) = this->b;
        b.segment(this->nC, this->n) = x - this->c;

        // solve MIQP
        const OptSolution sol = this->mi_opt(P, q, 0, A, b, settings, solution, warm_start_params);
        return !(sol.infeasible);
    }


    OptSolution HybZono::mi_opt(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                                const zono_float c, const Eigen::SparseMatrix<zono_float>& A,
                                const Eigen::Vector<zono_float, -1>& b,
                                const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                const WarmStartParams& warm_start_params) const
    {
        // QP data
        Eigen::Vector<zono_float, -1> xi_lb(this->nG);
        if (this->zero_one_form)
            xi_lb.setZero();
        else
            xi_lb.setConstant(-1);
        const Eigen::Vector<zono_float, -1> xi_ub = Eigen::Vector<zono_float, -1>::Ones(this->nG);

        const auto admm_data = std::make_shared<ADMM_data>(P, q, A, b, xi_lb, xi_ub, c, settings);

        // mixed integer data
        MI_data mi_data;
        mi_data.admm_data = admm_data;
        mi_data.idx_b = std::make_pair(this->nGc, this->nGb);
        mi_data.zero_one_form = this->zero_one_form;

        // build MI_ADMM_solver object
        MI_Solver mi_solver(mi_data);

        // warm start if applicable
        if (warm_start_params.z.size() == this->nG)
        {
            mi_solver.warmstart(warm_start_params.z, warm_start_params.u);
        }

        // solve optimization problem
        OptSolution sol = mi_solver.solve();

        if (solution != nullptr)
            *solution = std::make_shared<OptSolution>(sol);
        return sol;
    }

    std::vector<OptSolution> HybZono::mi_opt_multisol(const Eigen::SparseMatrix<zono_float>& P,
                                                      const Eigen::Vector<zono_float, -1>& q,
                                                      const zono_float c, const Eigen::SparseMatrix<zono_float>& A,
                                                      const Eigen::Vector<zono_float, -1>& b, int n_sols,
                                                      const OptSettings& settings,
                                                      std::shared_ptr<OptSolution>* solution) const
    {
        // ADMM data
        Eigen::Vector<zono_float, -1> xi_lb(this->nG);
        if (this->zero_one_form)
            xi_lb.setZero();
        else
            xi_lb.setConstant(-1);
        const Eigen::Vector<zono_float, -1> xi_ub = Eigen::Vector<zono_float, -1>::Ones(this->nG);

        const auto admm_data = std::make_shared<ADMM_data>(P, q, A, b, xi_lb, xi_ub, c, settings);

        // mixed integer data
        MI_data mi_data;
        mi_data.admm_data = admm_data;
        mi_data.idx_b = std::make_pair(this->nGc, this->nGb);
        mi_data.zero_one_form = this->zero_one_form;

        // build MI_ADMM_solver object
        MI_Solver mi_solver(mi_data);

        // solve optimization problem
        auto [fst, snd] = mi_solver.multi_solve(n_sols);
        if (solution != nullptr)
            *solution = std::make_shared<OptSolution>(snd);
        return fst;
    }


    void HybZono::remove_generators(Eigen::SparseMatrix<zono_float>& G, Eigen::SparseMatrix<zono_float>& A,
                                    const std::set<int>& idx_to_remove)
    {
        // declare triplets
        std::vector<Eigen::Triplet<zono_float>> triplets;

        // update G
        int delta_ind = 0;
        for (int k = 0; k < G.outerSize(); k++)
        {
            if (idx_to_remove.count(k))
            {
                ++delta_ind;
            }
            else
            {
                for (Eigen::SparseMatrix<zono_float>::InnerIterator it(G, k); it; ++it)
                {
                    triplets.emplace_back(static_cast<int>(it.row()), k - delta_ind, it.value());
                }
            }
        }
        G.resize(G.rows(), G.cols() - delta_ind);
        G.setFromTriplets(triplets.begin(), triplets.end());

        // update A
        triplets.clear();
        delta_ind = 0;
        for (int k = 0; k < A.outerSize(); k++)
        {
            if (idx_to_remove.count(k))
            {
                ++delta_ind;
            }
            else
            {
                for (Eigen::SparseMatrix<zono_float>::InnerIterator it(A, k); it; ++it)
                {
                    triplets.emplace_back(static_cast<int>(it.row()), k - delta_ind, it.value());
                }
            }
        }
        A.resize(A.rows(), A.cols() - delta_ind);
        A.setFromTriplets(triplets.begin(), triplets.end());
    }

    std::set<int> HybZono::find_unused_generators(const Eigen::SparseMatrix<zono_float>& G,
                                                  const Eigen::SparseMatrix<zono_float>& A)
    {
        std::set<int> idx_no_cons;
        for (int k = 0; k < A.outerSize(); k++)
        {
            bool is_unused = true;
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(A, k); it; ++it)
            {
                if (std::abs(it.value()) > zono_eps)
                {
                    is_unused = false;
                    break;
                }
            }

            if (is_unused)
            {
                idx_no_cons.insert(k);
            }
        }

        // check if any of idx_no_cons multiply only zeros
        std::set<int> idx_to_remove;
        for (int idx_no_con : idx_no_cons)
        {
            bool is_zero = true;
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(G, idx_no_con); it; ++it)
            {
                if (std::abs(it.value()) > zono_eps)
                {
                    is_zero = false;
                    break;
                }
            }

            if (is_zero)
            {
                idx_to_remove.insert(idx_no_con);
            }
        }

        return idx_to_remove;
    }

    void HybZono::make_G_A()
    {
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        get_triplets_offset<zono_float>(this->Gc, tripvec, 0, 0);
        get_triplets_offset<zono_float>(this->Gb, tripvec, 0, this->nGc);
        this->G.resize(this->n, this->nGc + this->nGb);
        this->G.setFromTriplets(tripvec.begin(), tripvec.end());

        tripvec.clear();
        get_triplets_offset<zono_float>(this->Ac, tripvec, 0, 0);
        get_triplets_offset<zono_float>(this->Ab, tripvec, 0, this->nGc);
        this->A.resize(this->nC, this->nGc + this->nGb);
        this->A.setFromTriplets(tripvec.begin(), tripvec.end());

        this->nG = this->nGc + this->nGb;
    }

    void HybZono::set_Ac_Ab_from_A()
    {
        std::vector<Eigen::Triplet<zono_float>> triplets_Ac, triplets_Ab;

        // iterate over A
        for (int k = 0; k < this->A.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->A, k); it; ++it)
            {
                if (it.col() < this->nGc)
                {
                    triplets_Ac.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
                }
                else
                {
                    triplets_Ab.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()) - this->nGc,
                                             it.value());
                }
            }
        }

        // set Ac, Ab
        this->Ac.resize(this->nC, this->nGc);
        this->Ac.setFromTriplets(triplets_Ac.begin(), triplets_Ac.end());
        this->Ab.resize(this->nC, this->nGb);
        this->Ab.setFromTriplets(triplets_Ab.begin(), triplets_Ab.end());
    }

    std::vector<Eigen::Vector<zono_float, -1>> HybZono::get_bin_leaves(const OptSettings& settings,
                                                                       std::shared_ptr<OptSolution>* solution,
                                                                       const int n_leaves) const
    {
        // optimize over P=I, q=0
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        P.setIdentity();
        const Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(this->nG);

        // solve
        std::vector<OptSolution> sol_vec = this->mi_opt_multisol(P, q, 0, this->A, this->b, n_leaves, settings,
                                                                 solution);

        // get leaves as conzonos
        std::vector<Eigen::Vector<zono_float, -1>> bin_leaves;
        for (auto& sol : sol_vec)
        {
            bin_leaves.emplace_back(sol.z.segment(this->nGc, this->nGb));
        }

        return bin_leaves;
    }

    std::vector<ConZono> HybZono::get_leaves(const bool remove_redundancy, const OptSettings& settings,
                                             std::shared_ptr<OptSolution>* solution, const int n_leaves,
                                             const int contractor_iter) const
    {
        // allocate all threads to branch and bound
        OptSettings settings_get_leaves = settings;
        settings_get_leaves.n_threads_bnb += settings.n_threads_admm_fp;
        settings_get_leaves.n_threads_admm_fp = 0;

        // get leaves as conzonos
        const std::vector<Eigen::Vector<zono_float, -1>> bin_leaves = this->get_bin_leaves(
            settings_get_leaves, solution, n_leaves);
        std::vector<ConZono> leaves;
        for (auto& xi_b : bin_leaves)
        {
            Eigen::Vector<zono_float, -1> cp = this->c + this->Gb * xi_b;
            Eigen::Vector<zono_float, -1> bp = this->b - this->Ab * xi_b;
            leaves.emplace_back(this->Gc, cp, this->Ac, bp, this->zero_one_form);
        }
        if (remove_redundancy)
        {
            for (auto& leaf : leaves)
            {
                leaf.remove_redundancy(contractor_iter);
            }
        }

        return leaves;
    }

    std::unique_ptr<ConZono> HybZono::convex_relaxation() const
    {
        if (this->is_empty_set())
        {
            return std::make_unique<EmptySet>(this->n);
        }
        else if (this->nG == 0)
        {
            return std::make_unique<Point>(this->c);
        }
        else if (this->nC == 0 && this->nGb == 0)
        {
            return std::make_unique<Zono>(this->G, this->c, this->zero_one_form);
        }
        else
        {
            return std::make_unique<ConZono>(this->G, this->c, this->A, this->b, this->zero_one_form);
        }
    }

    std::unique_ptr<HybZono> vrep_2_hybzono(const std::vector<Eigen::Matrix<zono_float, -1, -1>>& Vpolys,
                                            const bool expose_indicators)
    {
        // error handling
        if (Vpolys.empty())
        {
            throw std::invalid_argument("set_from_vrep: Vpolys must have at least one polytope.");
        }

        // dimensions
        const int n_polys = static_cast<int>(Vpolys.size());
        const int n_dims = static_cast<int>(Vpolys[0].cols());
        int n_verts; // declare

        // check if all polytopes have the same number of dimensions
        for (const auto& Vpoly : Vpolys)
        {
            if (Vpoly.cols() != n_dims)
            {
                throw std::invalid_argument("set_from_vrep: all polytopes must have the same number of dimensions.");
            }
        }

        // initialize V and M matrices as std::vectors
        // each entry is a row
        std::vector<Eigen::Matrix<zono_float, 1, -1>> V_vec, M_vec;
        Eigen::Matrix<zono_float, 1, -1> M_row(n_polys);

        // loop through each polytope
        for (int i = 0; i < n_polys; i++)
        {
            n_verts = static_cast<int>(Vpolys[i].rows());
            for (int j = 0; j < n_verts; j++)
            {
                // check if the vertex is already in V_vec
                auto vertex_equal = [&](const Eigen::Matrix<zono_float, 1, -1>& v) -> bool
                {
                    return (v - Vpolys[i].row(j)).norm() < zono_eps;
                };

                if (auto it_V = std::find_if(V_vec.begin(), V_vec.end(), vertex_equal); it_V == V_vec.end())
                {
                    V_vec.emplace_back(Vpolys[i].row(j));
                    M_row.setZero();
                    M_row(i) = 1;
                    M_vec.push_back(M_row);
                }
                else
                {
                    const int idx = static_cast<int>(std::distance(V_vec.begin(), it_V));
                    M_vec[idx](i) = 1;
                }
            }
        }

        const int nV = static_cast<int>(V_vec.size()); // number of unique vertices

        // convert to Eigen matrices
        Eigen::Matrix<zono_float, -1, -1> V(n_dims, nV);
        Eigen::Matrix<zono_float, -1, -1> M(nV, n_polys);
        for (int i = 0; i < nV; i++)
        {
            V.col(i) = V_vec[i];
            M.row(i) = M_vec[i];
        }

        // directly build hybzono in [0,1] form

        // declare
        std::vector<Eigen::Triplet<zono_float>> tripvec;

        // output dimension
        int n_out = n_dims;
        if (expose_indicators)
            n_out += static_cast<int>(n_polys);

        // Gc = [V, 0]
        Eigen::SparseMatrix<zono_float> Gc = V.sparseView();
        Gc.conservativeResize(n_out, 2 * nV);

        // Gb = [0]
        tripvec.clear();
        for (int i = n_dims; i < n_out; ++i)
        {
            tripvec.emplace_back(i, i - n_dims, one);
        }
        Eigen::SparseMatrix<zono_float> Gb(n_out, n_polys);
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        Gb.setFromSortedTriplets(tripvec.begin(), tripvec.end());
#else
        Gb.setFromTriplets(tripvec.begin(), tripvec.end());
#endif

        // c = 0
        Eigen::Vector<zono_float, -1> c(n_out);
        c.setZero();

        // Ac = [1^T, 0^T;
        //       0^T, 0^T;
        //       I, diag[sum(M, 2)]]
        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ac(2 + nV, 2 * nV);
        Eigen::SparseMatrix<zono_float> I_nv(nV, nV);
        I_nv.setIdentity();
        for (int i = 0; i < nV; i++)
        {
            tripvec.emplace_back(0, i, one);
        }
        get_triplets_offset<zono_float>(I_nv, tripvec, 2, 0);
        Eigen::Vector<zono_float, -1> sum_M = M.rowwise().sum();
        for (int i = 0; i < nV; i++)
        {
            tripvec.emplace_back(2 + i, nV + i, sum_M(i));
        }
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        // Ab = [0^T;
        //       1^T;
        //       -M]
        Eigen::SparseMatrix<zono_float> Ab(2 + nV, n_polys);
        tripvec.clear();
        for (int i = 0; i < n_polys; i++)
        {
            tripvec.emplace_back(1, i, one);
        }
        Eigen::SparseMatrix<zono_float> mM_sp = -M.sparseView();
        get_triplets_offset<zono_float>(mM_sp, tripvec, 2, 0);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());

        // b = [1;
        //      1;
        //      0]
        Eigen::Vector<zono_float, -1> b(2 + nV);
        b.setZero();
        b(0) = 1;
        b(1) = 1;

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, true, true);
    }

    zono_float HybZono::do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings& settings,
                                   std::shared_ptr<OptSolution>* solution, const WarmStartParams& warm_start_params)
    {
        // check dimensions
        if (this->n != d.size())
        {
            throw std::invalid_argument("Support: inconsistent dimensions.");
        }

        // if sharp, can solve as convex optimization problem
        if (this->sharp)
        {
            const auto Zc = this->convex_relaxation();
            return Zc->support(d, settings, solution);
        }

        // cost
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // solve MIQP
        const OptSolution sol = this->mi_opt(P, q, 0, this->A, this->b, settings, solution, warm_start_params);

        // check feasibility and return solution
        if (sol.infeasible) // Z is empty
            throw std::invalid_argument("Support: infeasible");
        else
            return d.dot(this->G * sol.z + this->c);
    }

    Box HybZono::do_bounding_box(const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                 const WarmStartParams& warm_start_params)
    {
        // if sharp, compute from convex relaxation
        if (this->sharp)
        {
            const auto Z_CR = this->convex_relaxation();
            return Z_CR->bounding_box(settings, solution);
        }

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
        zono_float s_neg, s_pos;

        // build QP for ADMM
        const Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // mixed-integer solution

        // get support in all box directions
        auto sol_in = std::make_shared<OptSolution>();
        for (int i = 0; i < this->n; i++)
        {
            // negative direction

            // update QP
            d.setZero();
            d(i) = -1;
            q = -this->G.transpose() * d;

            // solve
            OptSolution sol = this->mi_opt(P, q, 0, this->A, this->b, settings, &sol_in, warm_start_params);
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_neg = -d.dot(this->G * sol.z + this->c);

            // positive direction

            // update QP
            d.setZero();
            d(i) = 1;
            q = -this->G.transpose() * d;

            // solve
            sol = this->mi_opt(P, q, 0, this->A, this->b, settings, &sol_in, warm_start_params);
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_pos = d.dot(this->G * sol.z + this->c);

            // store bounds
            box[i] = Interval(s_neg, s_pos);
        }

        return box;
    }

    std::unique_ptr<HybZono> HybZono::do_complement(const zono_float delta_m, const bool remove_redundancy,
                                                    const OptSettings& settings,
                                                    std::shared_ptr<OptSolution>* solution, const int n_leaves,
                                                    const int contractor_iter)
    {
        // make sure set in [-1,1] form
        if (this->is_0_1_form()) this->convert_form();

        // need to get leaves and do complement for each leaf if Z is a hybzono
        auto leaves = this->get_leaves(remove_redundancy, settings, solution, n_leaves, contractor_iter);
        if (leaves.empty())
        {
            throw std::runtime_error("HybZono complement: set is empty.");
        }
        std::vector<std::unique_ptr<HybZono>> complements; // init
        for (auto& leaf : leaves)
        {
            complements.emplace_back(leaf.complement(delta_m));
        }
        std::unique_ptr<HybZono> Z_out;
        for (auto& comp : complements)
        {
            if (!Z_out)
            {
                Z_out.reset(comp->clone());
            }
            else
            {
                Z_out = intersection(*Z_out, *comp);
            }
        }
        return Z_out;
    }

    // type checking
    bool HybZono::is_point() const
    {
        const auto PointCast = dynamic_cast<const Point*>(this);
        return PointCast != nullptr;
    }

    bool HybZono::is_zono() const
    {
        const auto ZonoCast = dynamic_cast<const Zono*>(this);
        const auto PointCast = dynamic_cast<const Point*>(this);
        return (ZonoCast != nullptr) && (PointCast == nullptr);
    }

    bool HybZono::is_conzono() const
    {
        const auto ConZonoCast = dynamic_cast<const ConZono*>(this);
        const auto ZonoCast = dynamic_cast<const Zono*>(this);
        const auto EmptySetCast = dynamic_cast<const EmptySet*>(this);

        return (ConZonoCast != nullptr) && (ZonoCast == nullptr) && (EmptySetCast == nullptr);
    }

    bool HybZono::is_hybzono() const
    {
        const auto HybZonoCast = dynamic_cast<const HybZono*>(this);
        const auto ConZonoCast = dynamic_cast<const ConZono*>(this);

        return (HybZonoCast != nullptr) && (ConZonoCast == nullptr);
    }

    bool HybZono::is_empty_set() const
    {
        const auto EmptySetCast = dynamic_cast<const EmptySet*>(this);
        return EmptySetCast != nullptr;
    }

}
