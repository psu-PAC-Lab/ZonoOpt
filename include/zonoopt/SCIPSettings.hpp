#ifndef ZONOOPT_SCIP_SETTINGS_HPP_
#define ZONOOPT_SCIP_SETTINGS_HPP_

/**
 * @file SCIPSettings.hpp
 * @brief Settings for the dynamically-loaded SCIP solver backend.
 *
 * Hybrid design (same as GurobiSettings): a small set of commonly-used parameters
 * are exposed as typed std::optional fields; rare parameters can be set via the
 * bool/int/longint/real/char/str_params escape-hatch maps keyed by SCIP's
 * documented parameter name (e.g., "limits/time", "numerics/feastol").
 *
 * Any field/map entry left unset is not applied, so SCIP's own default remains.
 *
 * Full parameter list: https://www.scipopt.org/doc/html/PARAMETERS.php
 */

#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "SolverDataStructures.hpp"

namespace ZonoOpt
{
    namespace detail
    {
        bool scip_available();
        const std::string& scip_unavailable_reason();
    }

    /**
     * @brief Settings for the dynamically-loaded SCIP solver backend.
     *
     * Pass an instance to a ZonoOpt optimization method to route through SCIP.
     * If the SCIP shared library cannot be loaded at runtime, set_default_solver_settings
     * with a SCIPSettings throws std::runtime_error; per-call passes silently fall back
     * to the internal solver.
     *
     * Notes:
     *  - SCIP's objective is linear by definition. Quadratic objectives are modeled
     *    via an epigraph reformulation (auxiliary variable + quadratic constraint),
     *    transparent to the caller.
     */
    struct SCIPSettings : SolverSettings
    {
        // ---- Limits ----
        std::optional<double> TimeLimit;       ///< wall-clock time limit in seconds ("limits/time")
        std::optional<double> MemLimit;        ///< memory limit in MB ("limits/memory")
        std::optional<int>    SolutionLimit;   ///< stop after this many MIP solutions ("limits/solutions")
        std::optional<int>    NodeLimit;       ///< max B&B nodes ("limits/totalnodes" → applied to "limits/nodes")

        // ---- Tolerances ----
        std::optional<double> MIPGap;          ///< relative MIP optimality gap ("limits/gap")
        std::optional<double> MIPGapAbs;       ///< absolute MIP optimality gap ("limits/absgap")
        std::optional<double> FeasibilityTol;  ///< feasibility tolerance ("numerics/feastol")

        // ---- Behavior ----
        std::optional<int>    Threads;         ///< worker threads ("parallel/maxnthreads"); SCIP default is sequential unless built with parallel
        std::optional<int>    Seed;            ///< random seed shift ("randomization/randomseedshift")

        // ---- Display ----
        /// SCIP verbosity 0..5 ("display/verblevel"); 0 = silent. SCIP default is 4.
        std::optional<int>    VerbLevel;

        // ---- Escape-hatch maps for any SCIP parameter not exposed above ----
        // Keys are SCIP's documented parameter names (slash-separated paths),
        // e.g., real_params["numerics/dualfeastol"] = 1e-9.
        std::map<std::string, bool>         bool_params;
        std::map<std::string, int>          int_params;
        std::map<std::string, long long>    longint_params;
        std::map<std::string, double>       real_params;
        std::map<std::string, char>         char_params;
        std::map<std::string, std::string>  str_params;

        // ---- SolverSettings interface ------------------------------------------
        std::unique_ptr<SolverSettings> clone() const override
        {
            return std::make_unique<SCIPSettings>(*this);
        }

        void verify_available() const override
        {
            if (!detail::scip_available())
            {
                throw std::runtime_error("SCIPSettings: " + detail::scip_unavailable_reason());
            }
        }

        // ---- Debug print --------------------------------------------------------
        std::string print() const
        {
            std::stringstream ss;
            ss << "SCIPSettings (unset fields use SCIP defaults):\n";
            auto opt_int    = [&](const char* n, const std::optional<int>& v)    { if (v) ss << "  " << n << " = " << *v << "\n"; };
            auto opt_dbl    = [&](const char* n, const std::optional<double>& v) { if (v) ss << "  " << n << " = " << *v << "\n"; };

            opt_dbl("TimeLimit", TimeLimit);
            opt_dbl("MemLimit", MemLimit);
            opt_int("SolutionLimit", SolutionLimit);
            opt_int("NodeLimit", NodeLimit);
            opt_dbl("MIPGap", MIPGap);
            opt_dbl("MIPGapAbs", MIPGapAbs);
            opt_dbl("FeasibilityTol", FeasibilityTol);
            opt_int("Threads", Threads);
            opt_int("Seed", Seed);
            opt_int("VerbLevel", VerbLevel);

            for (const auto& kv : bool_params)    ss << "  [bool]    " << kv.first << " = " << (kv.second ? "true" : "false") << "\n";
            for (const auto& kv : int_params)     ss << "  [int]     " << kv.first << " = " << kv.second << "\n";
            for (const auto& kv : longint_params) ss << "  [longint] " << kv.first << " = " << kv.second << "\n";
            for (const auto& kv : real_params)    ss << "  [real]    " << kv.first << " = " << kv.second << "\n";
            for (const auto& kv : char_params)    ss << "  [char]    " << kv.first << " = '" << kv.second << "'\n";
            for (const auto& kv : str_params)     ss << "  [str]     " << kv.first << " = \"" << kv.second << "\"\n";
            return ss.str();
        }

        /**
         * @brief Construct a new SCIPSettings object
         * @throw std::runtime_error if the SCIP shared library cannot be dynamically loaded
         */
        SCIPSettings() { verify_available(); }
    };

    /**
     * @brief Solver-native solution metadata produced by the SCIP backend.
     *
     * When an OptSolution is produced by a SCIP solve, OptSolution::external_results
     * points to an instance of this class so the caller can read SCIP-specific status,
     * node counts, gap, and dual bound. Inspect via dynamic_cast / isinstance.
     */
    struct SCIPSolverResults : ExternalSolverResults
    {
        /// Raw SCIP_Status code (see SCIP's status enum; layout differs across versions).
        int status = 0;

        /// Total branch-and-bound nodes explored (SCIPgetNTotalNodes).
        long long node_count = 0;

        /// Relative gap at termination (SCIPgetGap). +infinity if no primal bound found.
        double mip_gap = 0.0;

        /// Best dual bound found (SCIPgetDualbound).
        double dual_bound = 0.0;

        /// Number of feasible solutions in SCIP's solution storage at termination (SCIPgetNSols).
        int n_sols_pool = 0;

        std::shared_ptr<ExternalSolverResults> clone() const override
        {
            return std::make_shared<SCIPSolverResults>(*this);
        }

        std::string print() const override
        {
            std::stringstream ss;
            ss << "SCIPSolverResults:\n";
            ss << "  status: " << status << "\n";
            ss << "  node_count: " << node_count << "\n";
            ss << "  mip_gap: " << mip_gap << "\n";
            ss << "  dual_bound: " << dual_bound << "\n";
            ss << "  n_sols_pool: " << n_sols_pool << "\n";
            return ss.str();
        }
    };
} // namespace ZonoOpt

#endif
