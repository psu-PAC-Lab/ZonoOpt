#pragma once

/**
 * @file ADMM.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Convex and mixed-integer ADMM implementations used within ZonoOpt.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <memory>
#include <atomic>
#include <mutex>
#include <random>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Box.hpp"
#include "CholeskyUtilities.hpp"
#include "SolverDataStructures.hpp"

/*
    Primary reference:
    Boyd, Stephen, et al.
    "Distributed optimization and statistical learning via the alternating direction method of multipliers."
    Foundations and Trends® in Machine learning 3.1 (2011): 1-122.
*/

namespace ZonoOpt
{
    namespace detail
    {
        /**
         * @brief Immutable, shareable core of an ADMM problem.
         *
         * Owns the heavy problem matrices and their LDLT factorizations. Every ADMM_data
         * produced by ADMM_data::clone() during a branch-and-bound / ADMM-FP solve points at
         * the same ADMM_prob, so the matrices and factors are stored exactly once per solve
         * instead of once per clone.
         *
         * Accessors are const and thread safe; the lazily-built members (A_rm, both LDLT
         * factorizations) are guarded by std::once_flag so concurrent B&B / ADMM-FP threads may
         * share one instance.
         */
        class ADMM_prob
        {
        public:
            ADMM_prob() = default;

            ADMM_prob(Eigen::SparseMatrix<zono_float> P, Eigen::SparseMatrix<zono_float> A,
                      Eigen::Vector<zono_float, -1> b, zono_float rho)
            {
                set(std::move(P), std::move(A), std::move(b), rho);
            }

            // always held through shared_ptr; never copied
            ADMM_prob(const ADMM_prob&) = delete;
            ADMM_prob& operator=(const ADMM_prob&) = delete;

            void set(Eigen::SparseMatrix<zono_float> P, Eigen::SparseMatrix<zono_float> A,
                     Eigen::Vector<zono_float, -1> b, zono_float rho);

            // eagerly built
            const Eigen::SparseMatrix<zono_float>& P() const { return m_P; }
            const Eigen::SparseMatrix<zono_float>& A() const { return m_A; }
            const Eigen::SparseMatrix<zono_float>& AT() const { return m_AT; }
            const Eigen::Vector<zono_float, -1>& b() const { return m_b; }
            int n_x() const { return m_n_x; }
            int n_cons() const { return m_n_cons; }
            zono_float rho() const { return m_rho; }

            // lazily built, thread safe
            const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm() const;
            const LDLT_solver& ldlt_M() const;
            const LDLT_solver& ldlt_AAT() const;

            // non-blocking readiness probes (verbosity/timing and tests only)
            bool M_ready() const { return m_M_ready.load(std::memory_order_acquire); }
            bool AAT_ready() const { return m_AAT_ready.load(std::memory_order_acquire); }
            bool A_rm_ready() const { return m_A_rm_ready.load(std::memory_order_acquire); }

        private:
            Eigen::SparseMatrix<zono_float> m_P, m_A, m_AT;
            Eigen::Vector<zono_float, -1> m_b;
            int m_n_x = 0, m_n_cons = 0;
            zono_float m_rho = 0;

            mutable Eigen::SparseMatrix<zono_float, Eigen::RowMajor> m_A_rm;
            mutable LDLT_solver m_ldlt_M, m_ldlt_AAT;
            mutable std::once_flag m_A_rm_once, m_M_once, m_AAT_once;
            mutable std::atomic<bool> m_M_ready{false}, m_AAT_ready{false}, m_A_rm_ready{false};

            void build_ldlt_M() const;
            void build_ldlt_AAT() const;
        };

        /**
         * @brief Per-instance ADMM state. Cheap to clone: shares the problem core (ADMM_prob).
         */
        struct ADMM_data
        {
            std::shared_ptr<const ADMM_prob> prob; // shared across all clones of one solve
            Eigen::Vector<zono_float, -1> q; // per-instance: mutated by ConZono::do_bounding_box
            Eigen::Vector<zono_float, 1> c;
            int n_x = 0, n_cons = 0; // mirrored from prob
            zono_float sqrt_n_x = 0;
            std::shared_ptr<Box> x_box;
            OptSettings settings;

            // constructor
            ADMM_data() = default;

            ADMM_data(Eigen::SparseMatrix<zono_float> P, Eigen::Vector<zono_float, -1> q,
                      Eigen::SparseMatrix<zono_float> A, Eigen::Vector<zono_float, -1> b,
                      const Eigen::Vector<zono_float, -1>& x_l, const Eigen::Vector<zono_float, -1>& x_u,
                      const zono_float c = 0, const OptSettings& settings = OptSettings())
            {
                set(std::move(P), std::move(q), std::move(A), std::move(b), x_l, x_u, c, settings);
            }

            // set method
            void set(Eigen::SparseMatrix<zono_float> P, Eigen::Vector<zono_float, -1> q,
                     Eigen::SparseMatrix<zono_float> A, Eigen::Vector<zono_float, -1> b,
                     const Eigen::Vector<zono_float, -1>& x_l, const Eigen::Vector<zono_float, -1>& x_u,
                     zono_float c = 0, const OptSettings& settings = OptSettings());

            // clone method: shares prob; deep-copies q, c, x_box, settings
            ADMM_data* clone() const;
        };

        // utilities
        class cycle_buffer
        {
        public:
            cycle_buffer(size_t N, zono_float eps_r);

            // returns false if value already contained in set
            bool insert(zono_float val);

            // get method
            const std::vector<zono_float>& get_vals() const;

        private:
            const size_t N;
            const zono_float eps_r;
            std::vector<zono_float> buffer;
        };

        /**
         * @brief Convex ADMM solver targeted at constrained zonotope optimization problems.
         *
         */
        class ADMM_solver
        {
        public:
            /**
             * @brief Construct a new admm solver object
             *
             * @param data
             */
            explicit ADMM_solver(const std::shared_ptr<ADMM_data>& data);

            /**
             * @brief Construct a new admm solver object
             *
             * @param other
             */
            ADMM_solver(const ADMM_solver& other);

            /**
             * @brief Destroy the admm solver object
             *
             */
            virtual ~ADMM_solver() = default;

            /**
             * @brief Warm-starts ADMM solver with primal and dual variables.
             *
             * @param x
             * @param u
             */
            virtual void warmstart(const Eigen::Vector<zono_float, -1>& x,
                                   const Eigen::Vector<zono_float, -1>& u);

            /**
             * @brief Optional pre-factorization of problem matrices.
             *
             * @throws std::runtime_error if factorization fails (typically because A is not full row rank).
             */
            virtual void factorize();

            /**
             * @brief Solves optimization problem using ADMM.
             *
             * @param stop
             * @return OptSolution
             * @throws std::invalid_argument if problem data dimensions are inconsistent.
             * @throws std::runtime_error if matrix factorization fails (e.g., A not full row rank).
             */
            OptSolution solve(std::atomic<bool>* stop);

            OptSolution solve();

        protected:
            // protected fields
            std::shared_ptr<ADMM_data> data;
            zono_float eps_prim = static_cast<zono_float>(1e-3), eps_dual = static_cast<zono_float>(1e-3);

            // startup method
            bool startup(Box& x_box, OptSolution& solution, const std::set<int>& contract_inds = std::set<int>());

            // core solve method
            virtual void solve_core(const Box& x_box, OptSolution& solution, std::atomic<bool>* stop);

            // warm start
            Eigen::Vector<zono_float, -1> x0, u0;

            // flags
            bool is_warmstarted = false;

            // check for infeasibility certificate
            bool is_infeasibility_certificate(const Eigen::Vector<zono_float, -1>& ek,
                                              const Eigen::Vector<zono_float, -1>& xk, const Box& x_box) const;

            bool check_problem_dimensions() const;
        };

        /**
         * @brief Mixed-integer ADMM solver with feasibility pump-inspired anti-cycling heuristic.
         *
         */
        class ADMM_FP_solver : public ADMM_solver
        {
        private:
            std::mt19937 rand_gen;

            void init_rand_gen();

        public:
            /**
             * @brief Construct a new ADMM solver feas pump object
             *
             * @param data
             */
            ADMM_FP_solver(const std::shared_ptr<ADMM_data>& data);

            /**
             * @brief Construct a new ADMM solver feas pump object
             *
             * @param other
             */
            ADMM_FP_solver(const ADMM_FP_solver& other) = default;

        protected:
            enum FP_Phase
            {
                Phase1,
                Phase2
            };

            // solve core method
            void solve_core(const Box& x_box, OptSolution& solution, std::atomic<bool>* stop) override;

            void ADMM_FP_core(const MI_Box& x_box, OptSolution& solution, std::atomic<bool>* stop, FP_Phase phase);
        };
    } // end namespace detail
} // end namespace ZonoOpt
