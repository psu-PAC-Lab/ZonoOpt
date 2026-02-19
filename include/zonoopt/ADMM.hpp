#ifndef ZONOOPT_ADMM_HPP_
#define ZONOOPT_ADMM_HPP_

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


#include <vector>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <memory>
#include <set>
#include <cmath>
#include <random>
#include <atomic>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "CholeskyUtilities.hpp"
#include "Intervals.hpp"
#include "SparseMatrixUtilities.hpp"
#include "SolverDataStructures.hpp"

/* 
    Primary reference: 
    Boyd, Stephen, et al. 
    "Distributed optimization and statistical learning via the alternating direction method of multipliers." 
    Foundations and Trends® in Machine learning 3.1 (2011): 1-122.
*/

namespace ZonoOpt
{
    /**
     * @brief Warm start parameters for optimization routines in ZonoOpt library.
     *
     */
    struct WarmStartParams
    {
        Eigen::Vector<zono_float, -1> z, u;
    };

    namespace detail
    {
        /**
         * @brief Data structure for ADMM solver.
         *
         */
        struct ADMM_data : std::enable_shared_from_this<ADMM_data>
        {
            Eigen::SparseMatrix<zono_float> P, A, AT;
            Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_rm;
            Eigen::Vector<zono_float, -1> q, b;
            Eigen::Vector<zono_float, 1> c;
            LDLT_data ldlt_data_M, ldlt_data_AAT;
            int n_x, n_cons;
            zono_float sqrt_n_x;
            std::shared_ptr<Box> x_box;
            OptSettings settings;

            // constructor
            ADMM_data() = default;

            ADMM_data(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                      const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                      const Eigen::Vector<zono_float, -1>& x_l, const Eigen::Vector<zono_float, -1>& x_u,
                      const zono_float c = 0, const OptSettings& settings = OptSettings())
            {
                set(P, q, A, b, x_l, x_u, c, settings);
            }

            // set method
            void set(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                     const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                     const Eigen::Vector<zono_float, -1>& x_l, const Eigen::Vector<zono_float, -1>& x_u,
                     const zono_float c = 0, const OptSettings& settings = OptSettings());

            // clone method
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
            explicit ADMM_solver(const ADMM_data& data);

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
             */
            virtual void factorize();

            /**
             * @brief Solves optimization problem using ADMM.
             *
             * @param stop
             * @return OptSolution
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

            // factor problem data
            void factorize_M() const;

            void factorize_AAT() const;

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
             * @param data
             */
            ADMM_FP_solver(const ADMM_data& data);

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

#endif
