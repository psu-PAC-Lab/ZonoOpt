#ifndef ZONOOPT_SOLVER_DATA_STRUCUTURES_HPP_
#define ZONOOPT_SOLVER_DATA_STRUCUTURES_HPP_

/**
 * @file SolverDataStructures.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Optimization settings and solution data structures for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <Eigen/Dense>
#include <memory>
#include <sstream>
#include <string>

namespace ZonoOpt
{
    /**
     * @brief Abstract base for all solver settings.
     *
     * The dynamic type of the SolverSettings instance passed to ZonoOpt's optimization
     * routines selects which solver backend is used:
     *   - OptSettings      → internal ADMM / branch-and-bound (default)
     *   - GurobiSettings   → dynamically-loaded Gurobi; silent fallback to OptSettings()
     *                        if Gurobi cannot be loaded at runtime
     *
     * Derived classes implement clone() so the program-wide default settings can hold
     * a polymorphic copy (see get_default_solver_settings / set_default_solver_settings).
     */
    struct SolverSettings
    {
        virtual ~SolverSettings() = default;

        /**
         * @brief Polymorphic copy. Must be overridden by every concrete subclass.
         */
        virtual std::unique_ptr<SolverSettings> clone() const = 0;

        /**
         * @brief Verify the solver backend selected by this settings type is usable.
         *
         * Called by set_default_solver_settings(); should throw std::runtime_error if
         * the backend cannot be initialized (e.g., GurobiSettings throws if the Gurobi
         * library cannot be dynamically loaded). Default implementation is a no-op,
         * appropriate for the internal solver which is always available.
         *
         * @throws std::runtime_error if the backend cannot be initialized (subclass implementations only; the default no-op does not throw).
         */
        virtual void verify_available() const {}

        /**
         * @brief Return name of the solver backend selected by this settings type, e.g., "ZonoOpt", "Gurobi", "SCIP".
         * 
         * @return std::string 
         */
        virtual std::string solver_name() const = 0;
    };

    /**
     * @brief Settings for the internal ZonoOpt ADMM / branch-and-bound solver.
     *
     */
    struct OptSettings : SolverSettings
    {
        // general settings

        /// display optimization progress
        bool verbose = false;

        /// print every verbose_interval iterations
        int verbosity_interval = 100;

        /// max time for optimization
        double t_max = std::numeric_limits<double>::max();

        // ADMM settings

        /// max convex admm iterations
        int k_max_admm = 5000;

        /// admm penalty parameter, higher prioritizes feasibility during iterations, lower prioritizes optimality
        zono_float rho = static_cast<zono_float>(10.0);

        /// dual convergence tolerance
        zono_float eps_dual = static_cast<zono_float>(1e-2);

        /// primal convergence tolerance
        zono_float eps_prim = static_cast<zono_float>(1e-3);

        /// check infeasibility every k_inf_check iterations
        int k_inf_check = 10;

        /// use infinity norm for convergence check (if false, scaled 2-norm is used)
        bool inf_norm_conv = true;

        /// flag to use interval contractor for constraint tightening / implication
        bool use_interval_contractor = true;

        /// number of interval contractor iterations
        int contractor_iter = 1;

        // mixed integer settings

        /// 0: best first, 1: best dive
        int search_mode = 0;

        /// flag to perform solution polishing
        bool polish = true;

        /// dual residual convergence tolerance during branch and bound and search
        zono_float eps_dual_search = static_cast<zono_float>(1e-1);

        /// primal residual convergence tolerance during branch and bound and search
        zono_float eps_prim_search = static_cast<zono_float>(1e-2);

        /// relative convergence tolerance
        zono_float eps_r = static_cast<zono_float>(1e-2);

        /// absolute convergence tolerance
        zono_float eps_a = static_cast<zono_float>(1e-1);

        /// max number of branch-and-bound iterations
        int k_max_bnb = static_cast<int>(1e5);

        /// max threads for branch and bound
        int n_threads_bnb = 4;

        /// max threads for ADMM-FP
        int n_threads_admm_fp = 3;

        /// enables single-threaded ADMM-FP solution, overrides n_threads_bnb, n_threads_admm_fp
        bool single_threaded_admm_fp = false;

        /// terminate if more than this many nodes are in branch and bound queue
        int max_nodes = static_cast<int>(1e5);

        /// when applying interval contractor in branch and bound, this is how deep to search the constraint tree for affected variables
        int contractor_tree_search_depth = 10;

        /// enable perturbations in ADMM-FP
        bool enable_perturb_admm_fp = true;

        /// enable restarts (significant perturbations) in ADMM-FP
        bool enable_restart_admm_fp = true;

        /// max ADMM iterations for ADMM-FP phase 1 (objective included)
        int k_max_admm_fp_ph1 = 10000;

        /// max ADMM iterations for ADMM-FP phase 2 (no objective)
        int k_max_admm_fp_ph2 = 90000;

        /// in ADMM-FP, this is the max size of the buffer that checks for cycles
        int cycle_detection_buffer_size = 20;

        /// relative tolerance for cycle detection, triggers perturbation
        zono_float eps_perturb = static_cast<zono_float>(1e-3);

        /// perform restart operation if primal residual does not improve over this many iterations in ADMM-FP
        int k_restart = 5000;

        /// enable rng seed for ADMM-FP
        bool enable_rng_seed = false;

        /// rng seed for ADMM-FP
        unsigned int rng_seed = 0;

        // polymorphic copy
        std::unique_ptr<SolverSettings> clone() const override
        {
            return std::make_unique<OptSettings>(*this);
        }

        // validity check
        /**
         * @brief Checks whether settings struct is valid
         *
         * @return validity boolean
         */
        bool settings_valid() const
        {
            const bool general_valid = t_max > 0 && verbosity_interval > 0;
            const bool admm_valid = (rho > 0 && k_max_admm > 0 &&
                eps_dual >= 0 && eps_prim >= 0 && k_inf_check >= 0 && contractor_iter > 0);
            const bool mi_valid = (eps_r >= 0 && eps_a >= 0 && eps_dual_search > 0 && eps_prim_search > 0 &&
                k_max_bnb > 0 && n_threads_bnb >= 0 && n_threads_admm_fp >= 0 && max_nodes > 0 &&
                contractor_tree_search_depth > 0 &&
                k_max_admm_fp_ph1 > 0 && cycle_detection_buffer_size > 0 &&
                (search_mode == 0 || search_mode == 1) &&
                eps_perturb > 0);

            return (general_valid && admm_valid && mi_valid);
        }

        /**
         * @brief displays settings as string
         *
         * @return string
         */
        std::string print() const
        {
            std::stringstream ss;
            ss << "OptSettings structure: " << std::endl;
            ss << "  verbose: " << (verbose ? "true" : "false") << std::endl;
            ss << "  verbosity_interval: " << verbosity_interval << std::endl;
            ss << "  t_max: " << t_max << std::endl;
            ss << "  k_max_admm: " << k_max_admm << std::endl;
            ss << "  rho: " << rho << std::endl;
            ss << "  eps_dual: " << eps_dual << std::endl;
            ss << "  eps_prim: " << eps_prim << std::endl;
            ss << "  k_inf_check: " << k_inf_check << std::endl;
            ss << "  inf_norm_conv: " << (inf_norm_conv ? "true" : "false") << std::endl;
            ss << "  use_interval_contractor: " << (use_interval_contractor ? "true" : "false") << std::endl;
            ss << "  contractor_iter: " << contractor_iter << std::endl;
            ss << "  search_mode: " << search_mode << std::endl;
            ss << "  polish: " << polish << std::endl;
            ss << "  eps_dual_search: " << eps_dual_search << std::endl;
            ss << "  eps_prim_search: " << eps_prim_search << std::endl;
            ss << "  eps_r: " << eps_r << std::endl;
            ss << "  eps_a: " << eps_a << std::endl;
            ss << "  k_max_bnb: " << k_max_bnb << std::endl;
            ss << "  n_threads_bnb: " << n_threads_bnb << std::endl;
            ss << "  n_threads_admm_fp: " << n_threads_admm_fp << std::endl;
            ss << "  single_threaded_admm_fp: " << (single_threaded_admm_fp ? "true" : "false") << std::endl;
            ss << "  max_nodes: " << max_nodes << std::endl;
            ss << "  contractor_tree_search_depth: " << contractor_tree_search_depth << std::endl;
            ss << "  enable_perturb_admm_fp: " << (enable_perturb_admm_fp ? "true" : "false") << std::endl;
            ss << "  k_max_admm_fp_ph1: " << k_max_admm_fp_ph1 << std::endl;
            ss << "  k_max_admm_fp_ph2: " << k_max_admm_fp_ph2 << std::endl;
            ss << "  cycle_detection_buffer_size: " << cycle_detection_buffer_size << std::endl;
            ss << "  eps_perturb: " << eps_perturb << std::endl;
            ss << "  k_restart: " << k_restart << std::endl;
            ss << "  enable_rng_seed: " << (enable_rng_seed ? "true" : "false") << std::endl;
            ss << "  rng_seed: " << rng_seed << std::endl;
            ss << "  enable_restart_admm_fp: " << (enable_restart_admm_fp ? "true" : "false") << std::endl;
            return ss.str();
        }

        std::string solver_name() const override
        {
            return "ZonoOpt";
        }

        /**
         * @brief OptSettings constructor
         * 
         */
        OptSettings() = default;
    };

    /**
     * @brief Abstract base for external-solver-specific solution metadata.
     *
     * When an OptSolution is produced by an external solver (Gurobi, SCIP, ...),
     * OptSolution::external_results points to a polymorphic ExternalSolverResults
     * subclass carrying the solver-native status code, node count, gap, etc.
     * The user can dynamic_cast to the concrete subclass to extract that info:
     *
     *     if (auto* gr = dynamic_cast<const GurobiSolverResults*>(sol.external_results.get())) {
     *         std::cout << "Gurobi MIP nodes: " << gr->node_count << "\n";
     *     }
     */
    struct ExternalSolverResults
    {
        virtual ~ExternalSolverResults() = default;

        /// Polymorphic copy; required so callers can deep-copy an OptSolution.
        virtual std::shared_ptr<ExternalSolverResults> clone() const = 0;

        /// Human-readable summary of the solver-specific fields.
        virtual std::string print() const { return "ExternalSolverResults"; }
    };

    /**
     * @brief Solution data structure for optimization routines in ZonoOpt library.
     *
     * The fields z, J, run_time, converged, and infeasible are always populated
     * regardless of which solver computed the solution. Fields below labeled
     * "admm-specific" carry meaningful values only when the internal solver was
     * used; when an external solver (Gurobi, SCIP, ...) was used,
     * `external_results` points to the solver-native solution metadata and the
     * admm-specific fields are left at their defaults.
     */
    struct OptSolution
    {
        // general — always populated

        /// solution vector
        Eigen::Vector<zono_float, -1> z;

        /// objective
        zono_float J = -std::numeric_limits<zono_float>::infinity();

        /// time to compute solution (seconds)
        double run_time = 0.0;

        /// true if optimization has converged (proved optimal for external solvers; satisfied tolerances for internal)
        bool converged = false;

        /// true if optimization problem is provably infeasible
        bool infeasible = false;

        // admm-specific — meaningful only when the internal solver was used.

        /// time to factorize matrices and run interval contractors
        double startup_time = 0.0;

        /// number of iterations
        int iter = 0;

        /// ADMM primal variable, approximately equal to z when converged
        Eigen::Vector<zono_float, -1> x;

        /// ADMM dual variable
        Eigen::Vector<zono_float, -1> u;

        /// primal residual, corresponds to feasibility
        zono_float primal_residual = std::numeric_limits<zono_float>::infinity();

        /// dual residual, corresponds to optimality
        zono_float dual_residual = std::numeric_limits<zono_float>::infinity();

        // external-solver-specific — populated only when an external solver was used.

        /// Polymorphic solver-native solution metadata. nullptr when the internal
        /// solver was used; otherwise points to a GurobiSolverResults, SCIPSolverResults,
        /// etc. Inspect via dynamic_cast (C++) or isinstance (Python).
        std::shared_ptr<ExternalSolverResults> external_results;

        /**
         * @brief displays solution as string
         *
         * @return string
         */
        std::string print() const
        {
            std::stringstream ss;
            ss << "OptSolution structure:" << std::endl;
            ss << "  z: vector of length " << z.size() << std::endl;
            ss << "  J: " << J << std::endl;
            ss << "  run_time: " << run_time << std::endl;
            ss << "  startup_time: " << startup_time << std::endl;
            ss << "  iter: " << iter << std::endl;
            ss << "  converged: " << (converged ? "true" : "false") << std::endl;
            ss << "  infeasible: " << (infeasible ? "true" : "false") << std::endl;
            ss << "  x: vector of length " << x.size() << std::endl;
            ss << "  u: vector of length " << u.size() << std::endl;
            ss << "  primal_residual: " << primal_residual << std::endl;
            ss << "  dual_residual: " << dual_residual << std::endl;
            if (external_results)
                ss << "  external_results: " << external_results->print() << std::endl;
            return ss.str();
        }
    };

    // ---- Program-wide default solver settings -------------------------------
    //
    // ZonoOpt's optimization methods take a SolverSettings argument whose default
    // is whatever set_default_solver_settings() last installed. This lets a user
    // switch the entire library to (e.g.) Gurobi by calling
    //     set_default_solver_settings(GurobiSettings{});
    // once at program startup, instead of passing GurobiSettings to every call.
    //
    // Initial default is OptSettings() (internal ADMM/branch-and-bound).
    //
    // Not thread-safe: set_default_solver_settings should be called from a single
    // thread at startup before parallel optimization begins. Existing const
    // SolverSettings& references obtained from get_default_solver_settings() are
    // invalidated by a subsequent set_default_solver_settings() call.

    namespace detail
    {
        inline std::unique_ptr<SolverSettings>& default_settings_storage()
        {
            static std::unique_ptr<SolverSettings> p = std::make_unique<OptSettings>();
            return p;
        }
    } // namespace detail

    /**
     * @brief Returns a reference to the current program-wide default solver settings.
     *
     * The returned reference is valid until the next call to set_default_solver_settings().
     */
    inline const SolverSettings& get_default_solver_settings()
    {
        return *detail::default_settings_storage();
    }

    /**
     * @brief Replaces the program-wide default solver settings with a polymorphic copy of `settings`.
     *
     * After this call, any ZonoOpt optimization method invoked without an explicit
     * settings argument will use a freshly resolved reference to the new default.
     *
     * @throws std::runtime_error if the backend selected by `settings` cannot be
     *         initialized (e.g., a GurobiSettings whose Gurobi library cannot be
     *         dynamically loaded). The existing default is left unchanged on throw.
     */
    inline void set_default_solver_settings(const SolverSettings& settings)
    {
        settings.verify_available();
        detail::default_settings_storage() = settings.clone();
    }
} // end namespace ZonoOpt

#endif
