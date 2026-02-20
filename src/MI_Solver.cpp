#include "ZonoOpt.hpp"

namespace ZonoOpt::detail
{
    MI_Solver::MI_Solver(const MI_data& data) : data(data), node_queue(comp)
    {
        // check settings validity
        if (!this->data.admm_data->settings.settings_valid())
        {
            throw std::invalid_argument("MI_solver setup: invalid settings.");
        }

        // modify multithreading settings if necessary
        if (!this->data.admm_data->settings.single_threaded_admm_fp)
        {
            if (std::thread::hardware_concurrency() < 2)
            {
                std::stringstream ss;
                ss << "Branch and bound requires at least 2 threads (1 for branch and bound, 1 for convergence monitoring), but only "
                    << std::thread::hardware_concurrency()
                    << " are available.";

                throw std::runtime_error(ss.str());
            }

            while (this->data.admm_data->settings.n_threads_bnb + this->data.admm_data->settings.n_threads_admm_fp >
                   static_cast<int>(std::thread::hardware_concurrency()) - 1)
            {
                if (this->data.admm_data->settings.n_threads_admm_fp > 0)
                {
                    this->data.admm_data->settings.n_threads_admm_fp--;
                }
                else
                {
                    this->data.admm_data->settings.n_threads_bnb--;
                }
            }
        }
    }

    OptSolution MI_Solver::solve()
    {
        this->multi_sol = false;
        auto sol = solver_core();
        return std::get<OptSolution>(sol);
    }

    std::pair<std::vector<OptSolution>, OptSolution> MI_Solver::multi_solve(const int max_sols)
    {
        this->multi_sol = true;
        auto sol = solver_core(max_sols);
        return std::get<std::pair<std::vector<OptSolution>, OptSolution>>(sol);
    }

    std::unique_ptr<Node, MI_Solver::NodeDeleter> MI_Solver::make_node(const std::shared_ptr<ADMM_data>& admm_data)
    {
        void* mem = pool.allocate(sizeof(Node), alignof(Node));
        auto node = new(mem) Node(admm_data);
        return {node, NodeDeleter(&pool)};
    }

    std::unique_ptr<Node, MI_Solver::NodeDeleter> MI_Solver::clone_node(const std::unique_ptr<Node, NodeDeleter>& other)
    {
        void* mem = pool.allocate(sizeof(Node), alignof(Node));
        Node* node = new(mem) Node(*other);
        return {node, NodeDeleter(&pool)};
    }

    std::variant<OptSolution, std::pair<std::vector<OptSolution>, OptSolution>> MI_Solver::solver_core(int max_sols)
    {
        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        double run_time;

        // set flags
        this->done = false;
        this->converged = false;

        // verbosity
        std::stringstream ss;
        if (this->data.admm_data->settings.verbose)
        {
            if (this->multi_sol)
            {
                ss << "Finding up to " << max_sols << " solutions to MIQP with " << this->data.admm_data->n_x <<
                    " variables and "
                    << this->data.admm_data->n_cons << " constraints using " << this->data.admm_data->settings.
                    n_threads_bnb << " branch-and-bound threads and "
                    << this->data.admm_data->settings.n_threads_admm_fp << " ADMM-FP threads.";
            }
            else if (this->data.admm_data->settings.single_threaded_admm_fp)
            {
                ss << "Searching for feasible solution to MIQP with " << this->data.admm_data->n_x << " variables and "
                    << this->data.admm_data->n_cons << " constraints using ADMM-FP (single-threaded).";
            }
            else
            {
                ss << "Solving MIQP problem with " << this->data.admm_data->n_x << " variables and "
                    << this->data.admm_data->n_cons << " constraints using " << this->data.admm_data->settings.
                    n_threads_bnb << " branch-and-bound threads and "
                    << this->data.admm_data->settings.n_threads_admm_fp << " ADMM-FP threads.";
            }
            print_str(ss);
        }

        auto return_infeasible_solution = [this, &start
            ]() -> std::variant<OptSolution, std::pair<std::vector<OptSolution>, OptSolution>>
        {
            OptSolution infeasible_solution;
            infeasible_solution.infeasible = true;
            infeasible_solution.run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<
                std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());
            infeasible_solution.startup_time = infeasible_solution.run_time;
            infeasible_solution.iter = 0;
            infeasible_solution.converged = false;
            infeasible_solution.J = std::numeric_limits<zono_float>::infinity();
            infeasible_solution.z = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.x = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.u = Eigen::Vector<zono_float, -1>::Zero(this->data.admm_data->n_x);
            infeasible_solution.primal_residual = std::numeric_limits<zono_float>::infinity();
            infeasible_solution.dual_residual = std::numeric_limits<zono_float>::infinity();

            if (this->multi_sol)
            {
                return std::make_pair(this->solutions.get(), infeasible_solution);
            }
            else
            {
                return infeasible_solution;
            }
        };

        // add root node
        this->bnb_data.reset(this->data.admm_data->clone()); // init
        this->bnb_data->settings.verbose = this->data.admm_data->settings.verbose &&
            this->data.admm_data->settings.single_threaded_admm_fp;
        this->bnb_data->settings.eps_dual = this->data.admm_data->settings.eps_dual_search;
        this->bnb_data->settings.eps_prim = this->data.admm_data->settings.eps_prim_search;
        std::unique_ptr<Node, NodeDeleter> root = this->make_node(this->bnb_data);
        if (!root->run_contractor()) // check for infeasibility during setup via interval contractor
        {
            return return_infeasible_solution();
        }

        // log startup time
        double startup_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        if (this->data.admm_data->settings.verbose)
        {
            ss << "Startup time = " << startup_time << " sec";
            print_str(ss);
        }

        // check if doing ADMM-FP warmstart
        const bool admm_fp_ws = this->data.admm_data->settings.single_threaded_admm_fp &&
            this->xi_ws.size() == this->data.admm_data->n_x;

        // solve root relaxation
        if (!admm_fp_ws)
        {
            root->solve();
            if (root->solution.infeasible)
            {
                return return_infeasible_solution();
            }
        }

        // make ADMM-FP thread

        // if contractor has collapsed any integer values, use that information
        this->admm_fp_data.reset(this->bnb_data->clone()); // copies over matrix factorization
        this->admm_fp_data->x_box = std::make_shared<MI_Box>(root->get_box().lower(), root->get_box().upper(),
                                                             this->data.idx_b, this->data.zero_one_form);
        const zono_float low = this->data.zero_one_form ? zero : -one;
        constexpr zono_float high = 1;
        for (int i = this->data.idx_b.first; i < this->data.idx_b.first + this->data.idx_b.second; i++)
        {
            bool low_cutoff = std::abs(root->get_box().lower()(i) - low) > zono_eps;
            bool high_cutoff = std::abs(root->get_box().upper()(i) - high) > zono_eps;
            if (low_cutoff && high_cutoff)
            {
                return return_infeasible_solution();
            }
            else if (low_cutoff)
            {
                (*this->admm_fp_data->x_box)[i] = Interval(high, high);
            }
            else if (high_cutoff)
            {
                (*this->admm_fp_data->x_box)[i] = Interval(low, low);
            }
        }

        // update ADMM-FP settings
        this->admm_fp_data->settings.use_interval_contractor = false;
        this->admm_fp_data->settings.eps_prim = this->data.admm_data->settings.eps_prim;
        this->admm_fp_data->settings.eps_dual = this->data.admm_data->settings.eps_dual;

        // make object
        std::unique_ptr<ADMM_FP_solver> pump = std::make_unique<ADMM_FP_solver>(this->admm_fp_data);

        // warm start with root relaxation solution
        if (admm_fp_ws)
        {
            const Eigen::Vector<zono_float, -1> u0 = this->u_ws.size() == this->data.admm_data->n_x
                                                         ? this->u_ws
                                                         : Eigen::Vector<zono_float, -1>::Zero(
                                                             this->data.admm_data->n_x);
            pump->warmstart(this->xi_ws, u0);
        }
        else
        {
            pump->warmstart(root->solution.z, root->solution.u);
        }

        // single-threaded case
        if (this->data.admm_data->settings.single_threaded_admm_fp)
        {
            run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());
            this->admm_fp_data->settings.t_max = this->data.admm_data->settings.t_max - run_time; // reset max time

            // run ADMM-FP
            OptSolution sol = pump->solve(nullptr);
            const int admm_fp_iter = sol.iter;

            // optionally polish
            if (sol.converged && this->data.admm_data->settings.polish)
            {
                const std::shared_ptr<ADMM_data> convex_node_data(this->bnb_data->clone());
                // copies over matrix factorization
                for (int i = this->data.idx_b.first; i < this->data.idx_b.first + this->data.idx_b.second; i++)
                {
                    (*convex_node_data->x_box)[i] = Interval(sol.z(i), sol.z(i));
                }
                convex_node_data->settings.eps_prim = this->data.admm_data->settings.eps_prim; // strengthen tolerances
                convex_node_data->settings.eps_dual = this->data.admm_data->settings.eps_dual; // strengthen tolerances
                convex_node_data->settings.use_interval_contractor = true; // use interval contractor

                run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start).count());
                convex_node_data->settings.t_max = this->data.admm_data->settings.t_max - run_time;

                ADMM_solver convex_node(convex_node_data);
                convex_node.warmstart(sol.z, sol.u);
                sol = convex_node.solve(nullptr);
            }

            // solution time
            sol.infeasible = false; // cannot prove infeasibility here
            sol.iter = admm_fp_iter; // log the admm-fp iterations
            sol.run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());

            return sol;
        }
        // start threads
        std::vector<std::thread> bnb_threads;
        for (int i = 0; i < this->data.admm_data->settings.n_threads_bnb; i++)
        {
            bnb_threads.emplace_back([this]() { worker_loop(); });
        }

        std::vector<std::thread> fp_threads;
        for (int i = 0; i < this->data.admm_data->settings.n_threads_admm_fp; i++)
        {
            auto admm_fp_node = std::make_unique<ADMM_FP_solver>(*pump);
            fp_threads.emplace_back([this, node=std::move(admm_fp_node)]() mutable { admm_fp_loop(std::move(node)); });
        }

        // push root to node queue
        this->push_node(std::move(root));

        // loop and check for exit conditions
        int print_iter = 0;
        if (this->data.admm_data->settings.verbose)
        {
            ss << std::endl << std::setw(13) << "Iter" << std::setw(13) << "Queue" << std::setw(13) <<
                "ADMM-FP Iter" << std::setw(13) <<
                "Threads" << std::setw(13) << "Time [s]" << std::setw(13) << "J_min" << std::setw(13) <<
                "J_max" << std::setw(13) << "Gap [%]" << std::setw(13) << "Feasible" << std::setw(13) <<
                "ADMM-FP sol" << std::endl;
            print_str(ss);
        }

        while (!this->done)
        {
            // check for timeout
            run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count());
            if (run_time > this->data.admm_data->settings.t_max)
            {
                this->done = true;
            }

            // check for max nodes
            int queue_size;
            {
                std::lock_guard<std::mutex> lock(pq_mtx);
                queue_size = static_cast<int>(this->node_queue.size());
            }
            if (queue_size > this->data.admm_data->settings.max_nodes)
            {
                this->done = true;
            }

            // check for max iterations
            if (this->iter >= this->data.admm_data->settings.k_max_bnb)
            {
                this->done = true;
            }

            // check for convergence

            // get lower bound / check if there are no nodes remaining
            zono_float J_min = -std::numeric_limits<zono_float>::infinity();
            {
                std::pair<zono_float, bool> J_min_threads_pair;
                std::lock_guard<std::mutex> lock(pq_mtx);
                J_min_threads_pair = this->J_threads.get_min(); // lower bound from active threads

                if (this->node_queue.empty())
                {
                    if (!J_min_threads_pair.second)
                    {
                        this->done = true; // no nodes remaining
                        this->converged = true;
                    }
                    else
                    {
                        J_min = J_min_threads_pair.first;
                    }
                }
                else
                {
                    J_min = std::min(this->node_queue.top()->solution.J, J_min_threads_pair.first);
                }
            }

            // check for convergence based on lower and upper bounds
            zono_float gap_percent = 0;
            if (!this->multi_sol)
            {
                const zono_float gap = std::abs(this->J_max - J_min);
                gap_percent = std::abs(this->J_max - J_min) / std::abs(this->J_max);
                if ((gap_percent < this->data.admm_data->settings.eps_r) || (gap < this->data.admm_data->settings.
                    eps_a))
                {
                    this->done = true;
                    this->converged = true;
                }
            }
            else // check based on number of solutions
            {
                if (this->solutions.size() >= static_cast<size_t>(max_sols))
                {
                    this->done = true;
                    this->converged = true;
                }
            }

            // verbosity
            if (this->data.admm_data->settings.verbose && (this->iter >= print_iter))
            {
                size_t n_threads = this->J_threads.size();
                ss << std::setw(13) << this->iter << std::setw(13) << queue_size << std::setw(13)
                    << this->iter_admm_fp << std::setw(13)
                    << n_threads << std::setw(13) << run_time << std::setw(13)
                    << J_min << std::setw(13) << this->J_max << std::setw(13)
                    << gap_percent * 100.0f << std::setw(13)
                    << (this->feasible ? "true" : "false") << std::setw(13)
                    << (this->feasible ? (this->admm_fp_incumbent ? "true" : "false") : "") << std::endl;
                print_str(ss);
                print_iter += this->data.admm_data->settings.verbosity_interval;
            }
        }

        // clean up
        pq_cv_bnb.notify_all(); // notify all threads to stop waiting
        pq_cv_admm_fp.notify_all();
        {
            std::lock_guard<std::mutex> lock(pq_mtx);
            this->node_queue.clear();
        }
        this->J_threads.clear();

        for (auto& thread : bnb_threads)
        {
            if (thread.joinable()) thread.join();
        }
        for (auto& thread : fp_threads)
        {
            if (thread.joinable()) thread.join();
        }

        this->node_queue.clear(); // nodes need to be freed before pool goes out of scope to avoid race condition

        // assemble solution
        OptSolution solution;
        solution.z = this->z.get();
        solution.J = this->J_max;
        solution.run_time = 1e-6 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count());
        solution.startup_time = startup_time;
        solution.iter = this->iter;
        solution.converged = this->converged;
        solution.infeasible = !this->feasible;

        solution.x = this->x.get();
        solution.u = this->u.get();
        solution.primal_residual = this->primal_residual;
        solution.dual_residual = this->dual_residual;

        if (this->multi_sol)
        {
            if (this->data.admm_data->settings.verbose)
            {
                ss << "Found " << this->solutions.size() << " solutions to MIQP problem.";
                print_str(ss);
            }

            return std::make_pair(solutions.get(), solution);
        }
        else
        {
            // verbosity
            if (this->data.admm_data->settings.verbose)
            {
                ss << "Converged = " << (this->converged ? "true" : "false") << ", Feasible = " <<
                    (this->feasible ? "true" : "false") << ", Iterations = " << this->iter << ", Solve time = " <<
                    solution.run_time << " sec, Objective = " << solution.J << ", Average QP iterations = " <<
                    static_cast<double>(this->qp_iter) / static_cast<double>(this->iter)
                    << ", Average solve time = " << this->total_run_time.get() / static_cast<double>(this->iter) <<
                    " sec, Average startup time = "
                    << this->total_startup_time.get() / static_cast<double>(this->iter) <<
                    " sec, Solution from ADMM-FP? "
                    << (this->admm_fp_incumbent ? "true" : "false") << std::endl;
                print_str(ss);
            }

            return solution;
        }
    }

    void MI_Solver::solve_and_branch(const std::unique_ptr<Node, NodeDeleter>& node)
    {
        // objective prior to solving
        const zono_float J_min_prior = node->solution.J;

        // solve node
        node->solve(&this->done);
        if (this->done) return;

        // cleanup function
        auto cleanup = [&, this]()
        {
            // remove J from J_threads vector
            this->J_threads.remove(J_min_prior);

            // increment nodes evaluated and logging info
            ++this->iter;
            this->qp_iter += node->solution.iter;
            this->total_run_time += node->solution.run_time;
            this->total_startup_time += node->solution.startup_time;
        };

        // return if infeasible, not converged, or no optimal solution exists in branch
        if (!(node->solution.infeasible || !node->solution.converged || (node->solution.J > this->J_max && !this->
            multi_sol)))
        {
            // check if node is integer feasible
            if (is_integer_feasible(node->solution.z.segment(this->data.idx_b.first, this->data.idx_b.second)))
            {
                // rerun with refined tolerance
                if (this->data.admm_data->settings.polish)
                {
                    node->update_convergence_tolerances(this->data.admm_data->settings.eps_prim,
                                                        this->data.admm_data->settings.eps_dual);
                    node->warmstart(node->solution.z, node->solution.u);
                    node->solve();
                }

                // make sure still integer feasible after refining
                if (!is_integer_feasible(node->solution.z.segment(this->data.idx_b.first, this->data.idx_b.second)))
                {
                    branch_most_frac(node);
                    cleanup(); // cleanup function
                    return;
                }

                // make sure return conditions are not met after refining
                if (node->solution.infeasible || !node->solution.converged || (node->solution.J > this->J_max && !this->
                    multi_sol))
                {
                    cleanup(); // cleanup function
                    return;
                }
                if (this->multi_sol) // store solution if doing multisol
                {
                    std::function<bool(const OptSolution&, const OptSolution&)> compare_eq = [this
                        ](const OptSolution& a, const OptSolution& b)
                    {
                        return this->check_bin_equal(a, b);
                    };
                    if (!this->solutions.contains(node->solution, compare_eq)) // new solution
                        this->solutions.push_back(node->solution);
                }
                if (node->solution.J < this->J_max - zono_eps) // check if node is better than current best
                {
                    // update incumbent
                    this->J_max = node->solution.J;
                    this->x.set(node->solution.x);
                    this->z.set(node->solution.z);
                    this->u.set(node->solution.u);
                    this->primal_residual = node->solution.primal_residual;
                    this->dual_residual = node->solution.dual_residual;
                    this->feasible = true;
                    this->admm_fp_incumbent = false;

                    // prune
                    if (!this->multi_sol) this->prune(node->solution.J);
                }
            }
            else
            {
                branch_most_frac(node);
            }
        }

        cleanup(); // cleanup function
    }

    void MI_Solver::admm_fp_solve(const std::unique_ptr<ADMM_FP_solver>& node)
    {
        // solve
        OptSolution sol = node->solve(&this->done);
        if (this->done)
        {
            ++this->iter_admm_fp;
            return;
        }

        // refine if converged and likely to be incumbent
        if (sol.converged && sol.J < this->J_max - zono_eps)
        {
            // polishing
            if (this->data.admm_data->settings.polish)
            {
                const std::shared_ptr<ADMM_data> convex_node_data(this->bnb_data->clone());
                // copies over matrix factorization
                for (int i = this->data.idx_b.first; i < this->data.idx_b.first + this->data.idx_b.second; i++)
                {
                    (*convex_node_data->x_box)[i] = Interval(sol.z(i), sol.z(i));
                }
                convex_node_data->settings.eps_prim = this->data.admm_data->settings.eps_prim; // strengthen tolerances
                convex_node_data->settings.eps_dual = this->data.admm_data->settings.eps_dual; // strengthen tolerances
                convex_node_data->settings.use_interval_contractor = true; // use interval contractor

                ADMM_solver convex_node(convex_node_data);
                convex_node.warmstart(sol.z, sol.u);
                sol = convex_node.solve(&this->done);
            }

            // make sure return conditions are not met after refining
            if (sol.infeasible || !sol.converged || (sol.J > this->J_max && !this->multi_sol))
            {
                ++this->iter_admm_fp;
                return;
            }

            // store solution if doing multisol
            if (this->multi_sol)
            {
                std::function<bool(const OptSolution&, const OptSolution&)> compare_eq = [this
                    ](const OptSolution& a, const OptSolution& b)
                {
                    return this->check_bin_equal(a, b);
                };
                if (!this->solutions.contains(sol, compare_eq)) // new solution
                    this->solutions.push_back(sol);
            }

            // new incumbent
            if (sol.J < this->J_max - zono_eps) // check if node is better than current best
            {
                // update incumbent
                this->J_max = sol.J;
                this->x.set(sol.x);
                this->z.set(sol.z);
                this->u.set(sol.u);
                this->primal_residual = sol.primal_residual;
                this->dual_residual = sol.dual_residual;
                this->feasible = true;
                this->admm_fp_incumbent = true;

                // prune
                if (!this->multi_sol) this->prune(sol.J);
            }
        }

        ++this->iter_admm_fp; // increment ADMM-FP iterations
    }

    bool MI_Solver::is_integer_feasible(const Eigen::Ref<const Eigen::Vector<zono_float, -1>> xb) const
    {
        const zono_float low = this->data.zero_one_form ? zero : -one;
        constexpr zono_float high = 1;

        for (int i = 0; i < xb.size(); i++)
        {
            if ((std::abs(xb(i) - high) > zono_eps) && (std::abs(xb(i) - low) > zono_eps))
                return false;
        }
        return true;
    }

    void MI_Solver::branch_most_frac(const std::unique_ptr<Node, NodeDeleter>& node)
    {
        // must be at least 1 binary variable
        if (this->data.idx_b.second <= 0)
            return;

        const zono_float low = this->data.zero_one_form ? zero : -one;
        constexpr zono_float high = 1;

        // round and find most fractional variable
        const Eigen::Array<zono_float, -1, 1> xb = node->solution.z.segment(
            this->data.idx_b.first, this->data.idx_b.second).array();
        Eigen::Array<zono_float, -1, 1> lower(xb.size());
        Eigen::Array<zono_float, -1, 1> upper(xb.size());
        lower.setConstant(low);
        upper.setConstant(high);
        const Eigen::Array<zono_float, -1, 1> d_l = (xb - lower).abs();
        const Eigen::Array<zono_float, -1, 1> d_u = (xb - upper).abs();
        const Eigen::Array<zono_float, -1, 1> d = d_l.min(d_u); // distance to rounded value
        int idx_most_frac;
        d.maxCoeff(&idx_most_frac); // index of most fractional variable
        idx_most_frac = this->data.idx_b.first + idx_most_frac; // convert to original index

        // branch on most fractional variable
        std::unique_ptr<Node, NodeDeleter> left = this->clone_node(node);
        std::unique_ptr<Node, NodeDeleter> right = this->clone_node(node);

        // branches
        const bool left_inf = !left->fix_bound(idx_most_frac, low);
        const bool right_inf = !right->fix_bound(idx_most_frac, high);

        // warm start
        left->warmstart(node->solution.z, node->solution.u);
        right->warmstart(node->solution.z, node->solution.u);

        switch (this->data.admm_data->settings.search_mode)
        {
        case (0):
            {
                // best first: push both nodes to queue
                if (!left_inf) this->push_node(std::move(left));
                if (!right_inf) this->push_node(std::move(right));
                break;
            }
        case (1):
            {
                // best dive: push worse nodes to queue, solve better node
                if (left_inf && right_inf) // both branches infeasible
                {
                    return; // nothing to do
                }
                else if (left_inf)
                {
                    this->solve_and_branch(right);
                }
                else if (right_inf)
                {
                    this->solve_and_branch(left);
                }
                else // both branches feasible
                {
                    if (left->get_box().width() > right->get_box().width()) // left is worse
                    {
                        this->push_node(std::move(left));
                        this->solve_and_branch(right);
                    }
                    else // right is worse
                    {
                        this->push_node(std::move(right));
                        this->solve_and_branch(left);
                    }
                }
                break;
            }
        default:
            {
                std::stringstream ss;
                ss << "MI_ADMM_solver: unknown search mode " << this->data.admm_data->settings.search_mode;
                throw std::runtime_error(ss.str());
            }
        }
    }

    void MI_Solver::worker_loop()
    {
        while (!this->done)
        {
            std::unique_ptr<Node, NodeDeleter> node = make_node(this->bnb_data);
            {
                std::unique_lock<std::mutex> lock(pq_mtx);
                pq_cv_bnb.wait(lock, [this]() { return this->done || !this->node_queue.empty(); });
                if (this->done) return;
                node = this->node_queue.pop_top();
                this->J_threads.add(node->solution.J);
                // add J to J_threads vector, need to do this before releasing lock
            }
            solve_and_branch(node);
        }
    }

    void MI_Solver::admm_fp_loop(std::unique_ptr<ADMM_FP_solver>&& node)
    {
        admm_fp_solve(node); // warm-started with root relaxation solution
        while (!this->done)
        {
            {
                std::unique_lock<std::mutex> lock(pq_mtx);
                pq_cv_admm_fp.wait(lock, [this]() { return this->done || !this->node_queue.empty(); });
                if (this->done) return;
                node->warmstart(this->node_queue.top()->solution.z, this->node_queue.top()->solution.u);
            }
            admm_fp_solve(node);
        }
    }

    void MI_Solver::push_node(std::unique_ptr<Node, NodeDeleter>&& node)
    {
        std::unique_lock<std::mutex> lock(pq_mtx);
        this->node_queue.push(std::move(node));
        pq_cv_bnb.notify_one();
        pq_cv_admm_fp.notify_one();
    }

    void MI_Solver::prune(const zono_float J_prune)
    {
        // create node with J = J_prune
        const std::unique_ptr<Node, NodeDeleter> n = this->make_node(this->bnb_data);
        n->solution.J = J_prune;
        {
            std::lock_guard<std::mutex> lock(pq_mtx);
            this->node_queue.prune(n);
        }
    }

    bool MI_Solver::check_bin_equal(const OptSolution& sol1, const OptSolution& sol2) const
    {
        return (sol1.z.segment(this->data.idx_b.first, this->data.idx_b.second) - sol2.z.segment(
                   this->data.idx_b.first, this->data.idx_b.second))
               .cwiseAbs().maxCoeff() < zono_eps;
    }
}
