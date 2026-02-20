#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    void ADMM_data::set(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                        const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                        const Eigen::Vector<zono_float, -1>& x_l, const Eigen::Vector<zono_float, -1>& x_u,
                        const zono_float c, const OptSettings& settings)
    {
        this->P = P;
        this->q = q;
        this->A = A;
        this->AT = A.transpose();
        this->A_rm = A;
        this->b = b;
        this->c(0) = c;

        this->n_x = static_cast<int>(P.rows());
        this->n_cons = static_cast<int>(A.rows());
        this->sqrt_n_x = std::sqrt(static_cast<zono_float>(this->n_x));

        this->x_box = std::make_shared<Box>(x_l, x_u);

        if (!settings.settings_valid()) throw std::invalid_argument("ADMM data: invalid settings.");
        this->settings = settings;
    }

    ADMM_data* ADMM_data::clone() const
    {
        const auto new_data = new ADMM_data(*this);
        new_data->x_box = std::make_shared<Box>(*this->x_box);
        return new_data;
    }

    ADMM_solver::ADMM_solver(const ADMM_data& data)
    {
        // store data and settings
        this->data = std::make_shared<ADMM_data>(data);
        this->eps_prim = data.settings.eps_prim;
        this->eps_dual = data.settings.eps_dual;

        // flags
        this->is_warmstarted = false;
    }

    ADMM_solver::ADMM_solver(const std::shared_ptr<ADMM_data>& data)
    {
        // store data and settings
        this->data = data;
        this->eps_prim = data->settings.eps_prim;
        this->eps_dual = data->settings.eps_dual;

        // flags
        this->is_warmstarted = false;
    }

    ADMM_solver::ADMM_solver(const ADMM_solver& other)
    {
        this->data = other.data;
        this->x0 = other.x0;
        this->u0 = other.u0;
        this->is_warmstarted = other.is_warmstarted;
        this->eps_dual = other.eps_dual;
        this->eps_prim = other.eps_prim;
    }

    void ADMM_solver::warmstart(const Eigen::Vector<zono_float, -1>& x,
                                const Eigen::Vector<zono_float, -1>& u)
    {
        // copy in warm start variables
        this->x0 = x;
        this->u0 = u;

        // set flag
        this->is_warmstarted = true;
    }

    void ADMM_solver::factorize()
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        double run_time;
        std::stringstream ss;
        if (!this->data->ldlt_data_M.factorized)
        {
            t0 = std::chrono::high_resolution_clock::now();
            this->factorize_M();
            if (this->data->settings.verbose)
            {
                run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - t0).count());
                ss << "M factorization time = " << run_time << " sec";
                print_str(ss);
            }
        }
        if (!this->data->ldlt_data_AAT.factorized)
        {
            t0 = std::chrono::high_resolution_clock::now();
            this->factorize_AAT();
            if (this->data->settings.verbose)
            {
                run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - t0).count());
                ss << "A*A^T factorization time = " << run_time << " sec";
                print_str(ss);
            }
        }
    }

    OptSolution ADMM_solver::solve(std::atomic<bool>* stop)
    {
        // declare solution
        OptSolution solution;

        // startup
        if (const bool contractor_feasible = this->startup(*this->data->x_box, solution); !contractor_feasible)
        {
            return solution;
        }

        // solve
        solve_core(*this->data->x_box, solution, stop);
        return solution;
    }

    OptSolution ADMM_solver::solve() { return this->solve(nullptr); }

    bool ADMM_solver::startup(Box& x_box, OptSolution& solution, const std::set<int>& contract_inds)
    {
        // start timer
        const auto start = std::chrono::high_resolution_clock::now();

        // reset verbosity
        std::stringstream ss;

        // check that problem data is consistent
        if (!this->check_problem_dimensions())
        {
            throw std::invalid_argument("ADMM solve: inconsistent problem data dimensions.");
        }
        if (this->data->settings.verbose)
        {
            ss << "Solving ADMM problem with " << this->data->n_x << " variables and " << this->data->n_cons <<
                " constraints.";
            print_str(ss);
        }

        // factorize if not already done
        this->factorize();

        // apply contractor
        bool contractor_feasible = true;
        if (this->data->settings.use_interval_contractor)
        {
            const auto t0 = std::chrono::high_resolution_clock::now();
            if (contract_inds.empty())
            {
                contractor_feasible = x_box.contract(this->data->A_rm, this->data->b,
                                                     this->data->settings.contractor_iter);
            }
            else
            {
                contractor_feasible = x_box.contract_subset(this->data->A_rm, this->data->b,
                                                            this->data->settings.contractor_iter,
                                                            this->data->A, contract_inds,
                                                            this->data->settings.contractor_tree_search_depth);
            }

            if (this->data->settings.verbose)
            {
                const double run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<
                    std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - t0).count());
                ss << "Interval contractor time = " << run_time << " sec";
                print_str(ss);
            }
        }

        // log startup time
        const double startup_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        solution.startup_time = startup_time;
        solution.run_time = startup_time;

        // early exit if contractor is infeasible
        if (!contractor_feasible)
        {
            if (this->data->settings.verbose)
            {
                ss << "Infeasibility detected via interval contractor";
                print_str(ss);
            }

            solution.infeasible = true;
            solution.iter = 0;
            solution.converged = false;
            solution.primal_residual = std::numeric_limits<zono_float>::infinity();
            solution.dual_residual = std::numeric_limits<zono_float>::infinity();
            solution.J = std::numeric_limits<zono_float>::infinity();
            solution.x = Eigen::Vector<zono_float, -1>::Zero(this->data->n_x);
            solution.z = solution.x; // 0
            solution.u = solution.x; // 0
        }
        return contractor_feasible;
    }

    void ADMM_solver::solve_core(const Box& x_box, OptSolution& solution, std::atomic<bool>* stop)
    {
        // start clock
        auto start = std::chrono::high_resolution_clock::now();
        std::stringstream ss;

        // initial values
        Eigen::Vector<zono_float, -1> xk, zk, uk, zkm1, rhs, x_nu;
        if (this->is_warmstarted)
        {
            xk = this->x0;
            uk = this->u0;
        }
        else
        {
            xk = 0.5 * (x_box.lower() + x_box.upper());
            uk = Eigen::Vector<zono_float, -1>::Zero(this->data->n_x);
        }
        zk = xk;
        rhs = Eigen::Vector<zono_float, -1>::Zero(this->data->n_x + this->data->n_cons);
        rhs.segment(this->data->n_x, this->data->n_cons) = this->data->b; // unchanging
        zkm1 = zk;

        // init residuals
        zono_float rp_k = std::numeric_limits<zono_float>::infinity(), rd_k = std::numeric_limits<
                       zono_float>::infinity();

        // init loop
        int k = 0;
        double run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        bool converged = false, infeasible = false;

        while ((k < this->data->settings.k_max_admm) && (run_time + solution.startup_time < this->data->settings.t_max)
            && !converged && !infeasible
            && !(stop && (*stop)))
        {
            // x update
            rhs.segment(0, this->data->n_x) = -this->data->q + this->data->settings.rho * (zk - uk);
            x_nu = solve_LDLT(this->data->ldlt_data_M, rhs);
            xk = x_nu.segment(0, this->data->n_x);

            // z update
            zk = xk + uk;
            x_box.project(zk);

            // u update
            uk += xk - zk;

            // check for infeasibility certificate
            if (k % this->data->settings.k_inf_check == 0)
            {
                infeasible = this->is_infeasibility_certificate(zk - xk, xk, x_box);
                if (this->data->settings.verbose && infeasible)
                {
                    ss << "Infeasibility certificate detected at iteration " << k;
                    print_str(ss);
                }
            }

            // check convergence
            if (this->data->settings.inf_norm_conv)
            {
                rp_k = (xk - zk).cwiseAbs().maxCoeff();
                rd_k = this->data->settings.rho * (zk - zkm1).cwiseAbs().maxCoeff();
                converged = (rp_k < this->eps_prim && rd_k < this->eps_dual);
            }
            else
            {
                rp_k = (xk - zk).norm();
                rd_k = this->data->settings.rho * (zk - zkm1).norm();
                converged = (rp_k < this->data->sqrt_n_x * this->eps_prim && rd_k < this->data->sqrt_n_x * this->
                    eps_dual);
            }

            // increment
            zkm1 = zk;
            ++k;

            // get time
            run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());

            // verbosity
            if (this->data->settings.verbose && (k % this->data->settings.verbosity_interval == 0))
            {
                ss << "k = " << k << ": primal residual = " << rp_k << ", dual residual = "
                    << rd_k << ", run time = " << run_time << " sec";
                print_str(ss);
            }
        }

        // verbosity
        if (this->data->settings.verbose)
        {
            if (converged)
            {
                ss << "ADMM converged in " << k << " iterations.";
                print_str(ss);
            }
            else if (infeasible)
            {
                ss << "ADMM detected infeasibility.";
                print_str(ss);
            }
            else
            {
                ss << "ADMM did not converge in " << k << " iterations.";
                print_str(ss);
            }
        }

        // reset flags
        this->is_warmstarted = false;

        // build output
        solution.x = xk;
        solution.z = zk;
        solution.u = uk;
        solution.J = (0.5 * zk.transpose() * this->data->P * zk + this->data->q.transpose() * zk + this->data->c)(0);
        solution.primal_residual = rp_k;
        solution.dual_residual = rd_k;
        solution.run_time = run_time + solution.startup_time;
        solution.iter = k;
        solution.converged = converged;
        solution.infeasible = infeasible;
    }

    void ADMM_solver::factorize_M() const
    {
        // system matrix
        Eigen::SparseMatrix<zono_float> M(this->data->n_x + this->data->n_cons, this->data->n_x + this->data->n_cons);

        Eigen::SparseMatrix<zono_float> I(this->data->n_x, this->data->n_x);
        I.setIdentity();
        Eigen::SparseMatrix<zono_float> Phi = this->data->P + this->data->settings.rho * I;

        std::vector<Eigen::Triplet<zono_float>> triplets;
        get_triplets_offset<zono_float>(Phi, triplets, 0, 0);
        get_triplets_offset<zono_float>(this->data->A, triplets, this->data->n_x, 0);
        get_triplets_offset<zono_float>(this->data->AT, triplets, 0, this->data->n_x);
        M.setFromTriplets(triplets.begin(), triplets.end());

        // factorize system matrix
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>> ldlt_solver_M;
        ldlt_solver_M.compute(M);
        if (ldlt_solver_M.info() != Eigen::Success)
            throw std::runtime_error("ADMM: factorization of problem data failed, most likely A is not full row rank");

        get_LDLT_data(ldlt_solver_M, this->data->ldlt_data_M);
    }

    void ADMM_solver::factorize_AAT() const
    {
        // factorize A*AT
        const Eigen::SparseMatrix<zono_float> AAT = this->data->A * this->data->AT;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>> ldlt_solver_AAT;
        ldlt_solver_AAT.compute(AAT);
        if (ldlt_solver_AAT.info() != Eigen::Success)
            throw std::runtime_error("ADMM: factorization of A*A^T failed, most likely A is not full row rank");
        get_LDLT_data(ldlt_solver_AAT, this->data->ldlt_data_AAT);
    }

    bool ADMM_solver::is_infeasibility_certificate(const Eigen::Vector<zono_float, -1>& ek,
                                                   const Eigen::Vector<zono_float, -1>& xk, const Box& x_box) const
    {
        // project ek onto row space of A (i.e. column space of AT)
        Eigen::Vector<zono_float, -1> A_e = this->data->A * ek;
        Eigen::Vector<zono_float, -1> AAT_inv_A_e = solve_LDLT(this->data->ldlt_data_AAT, A_e);
        Eigen::Vector<zono_float, -1> ek_proj = this->data->AT * AAT_inv_A_e;

        // check if this is an infeasibility certificate
        const zono_float e_x = ek_proj.dot(xk);
        const Interval e_box = x_box.dot(ek_proj);
        return !e_box.contains(e_x);
    }

    bool ADMM_solver::check_problem_dimensions() const
    {
        const bool prob_data_consistent = (this->data->P.rows() == this->data->n_x && this->data->P.cols() == this->data
            ->n_x &&
            this->data->q.size() == this->data->n_x && this->data->A.rows() == this->data->n_cons &&
            this->data->A.cols() == this->data->n_x && this->data->b.size() == this->data->n_cons &&
            this->data->x_box->size() == static_cast<size_t>(this->data->n_x));

        bool warm_start_consistent;
        if (this->is_warmstarted)
        {
            warm_start_consistent = (this->x0.size() == this->data->n_x && this->u0.size() == this->data->n_x &&
                this->data->x_box->size() == static_cast<size_t>(this->data->n_x));
        }
        else
            warm_start_consistent = true;

        return prob_data_consistent && warm_start_consistent;
    }

    cycle_buffer::cycle_buffer(const size_t N, const zono_float eps_r) : N(N), eps_r(eps_r)
    {
        buffer.reserve(N);
    }

    bool cycle_buffer::insert(const zono_float val)
    {
        for (const zono_float b : buffer)
        {
            if (std::abs(b - val) < eps_r * std::abs(val))
            {
                buffer.clear();
                return false;
            }
        }

        if (buffer.size() == N)
        {
            buffer.erase(buffer.begin());
        }
        buffer.push_back(val);

        return true;
    }

    const std::vector<zono_float>& cycle_buffer::get_vals() const
    {
        return buffer;
    }

    int perturb_binaries(const MI_Box& x_box, const Eigen::Ref<const Eigen::Vector<zono_float, -1>> xk,
                         Eigen::Ref<Eigen::Vector<zono_float, -1>> zk, const bool restart, std::mt19937& rand_gen)
    {
        // random flipping based on how fractional the variable is
        constexpr zono_float mp3 = static_cast<zono_float>(-0.3);
        constexpr zono_float p7 = static_cast<zono_float>(0.7);
        const zono_float dist_min = restart ? mp3 : zero;
        const zono_float dist_max = restart ? p7 : one;
        std::uniform_real_distribution<zono_float> distr(dist_min, dist_max);

        int bins_flipped = 0;
        for (int i = x_box.binary_indices().first; i < x_box.binary_indices().first + x_box.binary_indices().second; ++
             i)
        {
            // get fraction
            const zono_float w = x_box.upper()(i) - x_box.lower()(i);
            if (w < zono_eps) continue;

            // otherwise compute fraction
            const zono_float frac = std::abs(xk(i) - zk(i)) / w;

            // flip with probability proportional to fraction
            if (restart ? frac + std::max(distr(rand_gen), zero) > p5 : distr(rand_gen) < frac)
            {
                if (zk(i) < (x_box.lower()(i) + x_box.upper()(i)) * p5)
                {
                    zk(i) = x_box.upper()(i);
                }
                else
                {
                    zk(i) = x_box.lower()(i);
                }
                ++bins_flipped;
            }
        }

        return bins_flipped;
    }

    void ADMM_FP_solver::init_rand_gen()
    {
        std::random_device rd;
        const unsigned int seed = this->data->settings.enable_rng_seed ? this->data->settings.rng_seed : rd();
        this->rand_gen = std::mt19937(seed);
    }

    ADMM_FP_solver::ADMM_FP_solver(const std::shared_ptr<ADMM_data>& data)
        : ADMM_solver(data)
    {
        this->init_rand_gen();
    }

    ADMM_FP_solver::ADMM_FP_solver(const ADMM_data& data)
        : ADMM_solver(data)
    {
        this->init_rand_gen();
    }

    void ADMM_FP_solver::solve_core(const Box& x_box, OptSolution& solution, std::atomic<bool>* stop)
    {
        // cast to MI_Box
        const auto* mi_box = dynamic_cast<const MI_Box*>(&x_box);
        if (mi_box == nullptr)
        {
            ADMM_solver::solve_core(x_box, solution, stop);
            return; // use convex solver
        }

        // call feasibility pump core method
        ADMM_FP_core(*mi_box, solution, stop, Phase1);

        // run alternating projection method if not converged
        if (!solution.converged)
        {
            if (this->data->settings.verbose)
            {
                std::stringstream ss;
                ss << "Beginning Phase 2 of ADMM-FP";
                print_str(ss);
            }
            warmstart(solution.z, solution.u);
            ADMM_FP_core(*mi_box, solution, stop, Phase2);
        }
    }

    void ADMM_FP_solver::ADMM_FP_core(const MI_Box& x_box, OptSolution& solution, std::atomic<bool>* stop,
                                      FP_Phase phase)
    {
        // start clock
        const auto start = std::chrono::high_resolution_clock::now();
        std::stringstream ss;

        // initial values
        Eigen::Vector<zono_float, -1> xk, zk, zkm1, uk, rhs, x_nu;
        if (this->is_warmstarted)
        {
            xk = this->x0;
            uk = this->u0;
        }
        else
        {
            xk = 0.5 * (x_box.lower() + x_box.upper());
            uk = Eigen::Vector<zono_float, -1>::Zero(this->data->n_x);
        }
        zk = xk;
        rhs = Eigen::Vector<zono_float, -1>::Zero(this->data->n_x + this->data->n_cons);
        rhs.segment(this->data->n_x, this->data->n_cons) = this->data->b; // unchanging
        zkm1 = zk;

        // init residuals
        zono_float rp_k = std::numeric_limits<zono_float>::infinity(), rp_min = std::numeric_limits<
                       zono_float>::infinity(), rd_k = std::numeric_limits<zono_float>::infinity();
        cycle_buffer buf(this->data->settings.cycle_detection_buffer_size, this->data->settings.eps_perturb);

        // incumbent solution {residual, x, z, u}
        auto incumbent = std::make_tuple(rp_k, xk, zk, uk);

        // init loop
        int k = 0, k_restart = 0, restarts = 0;
        double run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        bool converged = false;
        int bins_flipped = 0;
        const int max_iter = phase == Phase1
                                 ? this->data->settings.k_max_admm_fp_ph1
                                 : this->data->settings.k_max_admm_fp_ph2;

        while ((k < max_iter) && (run_time + solution.run_time < this->data->settings.t_max) &&
            !converged && !(stop && (*stop)))
        {
            // x update
            switch (phase)
            {
            case Phase1:
                rhs.segment(0, this->data->n_x) = -this->data->q + this->data->settings.rho * (zk - uk);
                x_nu = solve_LDLT(this->data->ldlt_data_M, rhs);
                xk = x_nu.segment(0, this->data->n_x);
                break;
            case Phase2:
                xk = zk - uk;
                affine_set_projection(xk, this->data->A, this->data->AT, this->data->b, this->data->ldlt_data_AAT);
                break;
            default:
                throw std::invalid_argument("Unknown phase type");
            }

            // z update
            zk = xk + uk;
            x_box.project(zk);

            // u update
            uk += xk - zk;

            // compute residual and check for convergence
            if (this->data->settings.inf_norm_conv)
            {
                rp_k = (xk - zk).cwiseAbs().maxCoeff();
                converged = rp_k < this->eps_prim;
                rd_k = this->data->settings.rho * (zk - zkm1).cwiseAbs().maxCoeff();
                // only for populating output struct
            }
            else
            {
                rp_k = (xk - zk).norm();
                converged = rp_k < this->data->sqrt_n_x * this->eps_prim;
                rd_k = this->data->settings.rho * (zk - zkm1).norm(); // only for populating output struct
            }

            // update incumbent if residual improved
            if (rp_k < std::get<0>(incumbent) - zono_eps)
            {
                std::get<0>(incumbent) = rp_k;
                std::get<1>(incumbent) = xk;
                std::get<2>(incumbent) = zk;
                std::get<3>(incumbent) = uk;
            }

            // add primal residual to buffer and perturb binaries if cycle detected
            if (this->data->settings.enable_perturb_admm_fp && !buf.insert(rp_k))
            {
                bins_flipped += perturb_binaries(x_box, xk, zk, false, this->rand_gen);
            }

            // restart conditions
            if (this->data->settings.enable_restart_admm_fp)
            {
                if (rp_k < rp_min - zono_eps)
                {
                    rp_min = rp_k;
                    k_restart = 0;
                }
                else
                {
                    ++k_restart;
                }
                if (k_restart >= this->data->settings.k_restart)
                {
                    bins_flipped += perturb_binaries(x_box, xk, zk, true, this->rand_gen);
                    k_restart = 0;
                    rp_min = rp_k; // updated rp_k
                    ++restarts;
                }
            }

            // increment
            ++k;
            zkm1 = zk;

            // get time
            run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());

            // verbosity
            if (this->data->settings.verbose && (k % this->data->settings.verbosity_interval == 0))
            {
                ss << "k = " << k << ": primal residual = " << rp_k << ", run time = " << run_time << " sec, " <<
                    bins_flipped << " binaries perturbed, " << restarts << " restarts";
                print_str(ss);
                bins_flipped = 0;
            }
        }

        // verbosity
        if (this->data->settings.verbose)
        {
            if (converged)
            {
                ss << "ADMM-FP converged in " << k << " iterations.";
                print_str(ss);
            }
            else
            {
                ss << "ADMM-FP did not converge in " << k << " iterations.";
                print_str(ss);
            }
        }

        // reset flags
        this->is_warmstarted = false;

        // retrieve incumbent
        xk = std::get<1>(incumbent);
        zk = std::get<2>(incumbent);
        uk = std::get<3>(incumbent);

        // build output
        solution.x = xk;
        solution.z = zk;
        solution.u = uk;
        solution.J = (0.5 * zk.transpose() * this->data->P * zk + this->data->q.transpose() * zk + this->data->c)(0);
        solution.primal_residual = rp_k;
        solution.dual_residual = rd_k; // for logging only, not used for convergence check
        solution.run_time += run_time;
        solution.iter += k;
        solution.converged = converged;
        solution.infeasible = false; // infeasibility is not checked in feasibility pump
    }
}
