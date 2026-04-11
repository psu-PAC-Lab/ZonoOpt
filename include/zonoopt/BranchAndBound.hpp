#ifndef ZONOOPT_BRANCH_AND_BOUND_
#define ZONOOPT_BRANCH_AND_BOUND_

/**
 * @file BranchAndBound.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Internal branch-and-bound routines for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <sstream>
#include <iomanip>
#include <variant>
#include <memory_resource>
#include <cmath>
#include <random>

#include "BnbDataStructures.hpp"
#include "SolverDataStructures.hpp"
#include "ADMM.hpp"

namespace ZonoOpt::detail
{
    class BranchAndBound
    {
    public:
        explicit BranchAndBound(const MI_data& data);

        // solve
        OptSolution solve();

        // branch and bound where all possible solutions are returned
        std::pair<std::vector<OptSolution>, OptSolution> multi_solve(int max_sols = std::numeric_limits<int>::max());

        // warmstart - only used for ADMM-FP
        void warmstart(const Eigen::Vector<zono_float, -1>& xi_ws, const Eigen::Vector<zono_float, -1>& u_ws)
        {
            this->xi_ws = xi_ws;
            this->u_ws = u_ws;
        }

    private:
        struct NodeDeleter
        {
            std::pmr::synchronized_pool_resource* pool_ptr;

            explicit NodeDeleter(std::pmr::synchronized_pool_resource* pool_ptr) : pool_ptr(pool_ptr)
            {
            }

            void operator()(Node* node) const
            {
                if (node)
                {
                    node->~Node();
                    pool_ptr->deallocate(node, sizeof(Node), alignof(Node));
                }
            }
        };

        struct NodeCompare
        {
            bool operator()(const std::unique_ptr<Node, NodeDeleter>& n1,
                            const std::unique_ptr<Node, NodeDeleter>& n2) const
            {
                if (n2->is_priority() && !n1->is_priority())
                    return true;
                else if (!n2->is_priority() && n1->is_priority())
                    return false;
                else
                    return n1->solution.J > n2->solution.J;
            }
        };

        struct JThreadCompare
        {
            bool operator()(const std::pair<int, zono_float>& v1,
                            const std::pair<int, zono_float>& v2) const
            {
                if (v1.second != v2.second) return v1.second < v2.second;
                return v1.first < v2.first;
            }
        };

        template <typename T, typename Comp=std::less<T>>
        struct ThreadGuard
        {
            ThreadSafeSet<T, Comp>& thread_tags;
            T tag;
            bool specified = false;

            explicit ThreadGuard(ThreadSafeSet<T, Comp>& thread_tags): thread_tags(thread_tags)
            {
            }

            void specify_tag(const T& tag)
            {
                if (!specified)
                {
                    this->tag = tag;
                    this->specified = true;
                    this->thread_tags.add(tag);
                }
            }

            ~ThreadGuard()
            {
                if (specified)
                    this->thread_tags.remove(tag);
            }
        };

        const MI_data data;
        std::pmr::synchronized_pool_resource pool;
        NodeCompare comp;

        PriorityQueuePrunable<std::unique_ptr<Node, NodeDeleter>, NodeCompare> node_queue; // priority queue for nodes
        mutable std::mutex pq_mtx;
        mutable std::mutex incumbent_mtx; // guards atomic check-and-update of incumbent (J_max, z, x, u, etc.)
        std::condition_variable pq_cv_bnb, pq_cv_admm_fp;
        // condition variables for branch-and-bound and ADMM-FP threads

        bool multi_sol = false;
        std::shared_ptr<ADMM_data> bnb_data, admm_fp_data; // data for branch-and-bound threads / ADMM-FP threads

        std::atomic<bool> converged = false;
        std::atomic<bool> done = false;
        std::atomic<bool> feasible = false; // feasible solution found
        std::atomic<bool> admm_fp_incumbent = false; // incumbent is from ADMM FP
        std::atomic<long int> qp_iter = 0; // number of QP iterations
        std::atomic<int> iter = 0; // number of iterations
        std::atomic<int> iter_admm_fp = 0; // number of feasibility pump iterations
        std::atomic<zono_float> J_max = std::numeric_limits<zono_float>::infinity(); // upper bound
        ThreadSafeAccess<Eigen::Vector<zono_float, -1>> z, x, u; // solution vector
        std::atomic<zono_float> primal_residual = std::numeric_limits<zono_float>::infinity();
        std::atomic<zono_float> dual_residual = std::numeric_limits<zono_float>::infinity();
        ThreadSafeIncrementable<double> total_startup_time{0.0};
        ThreadSafeIncrementable<double> total_run_time{0.0};
        ThreadSafeSet<std::pair<int, zono_float>, JThreadCompare> J_threads; // threads for J values
        ThreadSafeVector<OptSolution> solutions; // solutions found
        std::uniform_int_distribution<int> uniform_dist{0, std::numeric_limits<int>::max()};
        std::mt19937 rng{0};

        // warmstart variables
        Eigen::Vector<zono_float, -1> xi_ws, u_ws;

        // allocate nodes
        std::unique_ptr<Node, NodeDeleter> make_node(const std::shared_ptr<ADMM_data>& admm_data);

        std::unique_ptr<Node, NodeDeleter> clone_node(const std::unique_ptr<Node, NodeDeleter>& other);

        // solver core
        std::variant<OptSolution, std::pair<std::vector<OptSolution>, OptSolution>> solver_core(
            int max_sols = std::numeric_limits<int>::max());

        // solve node and branch
        void solve_and_branch(const std::unique_ptr<Node, NodeDeleter>& node);

        void admm_fp_solve(const std::unique_ptr<ADMM_FP_solver>& node);

        // check if integer feasible, xb is vector of relaxed binary variables
        bool is_integer_feasible(const Eigen::Ref<const Eigen::Vector<zono_float, -1>> xb) const;

        // most fractional branching
        void branch_most_frac(const std::unique_ptr<Node, NodeDeleter>& node);

        // loop for multithreading
        void worker_loop();

        void admm_fp_loop(std::unique_ptr<ADMM_FP_solver>&& node);

        // push node to queue
        void push_node(std::unique_ptr<Node, NodeDeleter>&& node);

        // prune
        void prune(zono_float J_prune);

        // check if 2 solutions correspond to the same binaries
        bool check_bin_equal(const OptSolution& sol1, const OptSolution& sol2) const;
    };
} // namespace ZonoOpt::detail


#endif
