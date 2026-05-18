#ifndef ZONOOPT_GUROBI_SETTINGS_HPP_
#define ZONOOPT_GUROBI_SETTINGS_HPP_

/**
 * @file GurobiSettings.hpp
 * @brief Settings for the dynamically-loaded Gurobi solver backend.
 *
 * Hybrid design: commonly-used parameters are exposed as typed std::optional fields
 * (great IDE / autocompletion / docs), and rare parameters can be set via the
 * int_params / dbl_params / str_params escape-hatch maps keyed by Gurobi's
 * documented parameter name.
 *
 * Any field/map entry left unset (nullopt / absent) is not applied to the model,
 * leaving Gurobi's own default in place.
 *
 * Full parameter list and meanings: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html
 */

#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "SolverDataStructures.hpp"

namespace ZonoOpt
{
    // Forward declarations so GurobiSettings can check availability inline without
    // pulling in GurobiSolver.hpp (which depends on this file).
    namespace detail
    {
        bool gurobi_available();
        const std::string& gurobi_unavailable_reason();
    }

    /**
     * @brief Settings for the dynamically-loaded Gurobi solver backend.
     *
     * Pass an instance to a ZonoOpt optimization method to route through Gurobi.
     * If the Gurobi shared library cannot be loaded at runtime, the library silently
     * falls back to the internal solver with default OptSettings.
     */
    struct GurobiSettings : SolverSettings
    {
        // ---- Termination ----
        std::optional<double> TimeLimit;      ///< wall-clock time limit in seconds
        std::optional<double> WorkLimit;      ///< deterministic work limit
        std::optional<double> MemLimit;       ///< memory limit in GB
        std::optional<int>    SolutionLimit;  ///< stop after this many MIP solutions

        // ---- Tolerances ----
        std::optional<double> MIPGap;         ///< relative MIP optimality gap
        std::optional<double> MIPGapAbs;      ///< absolute MIP optimality gap
        std::optional<double> FeasibilityTol; ///< constraint feasibility tolerance
        std::optional<double> OptimalityTol;  ///< dual feasibility tolerance
        std::optional<double> IntFeasTol;     ///< integer feasibility tolerance

        // ---- Algorithm selection / behavior ----
        std::optional<int>    Method;         ///< root-node algorithm: -1 auto, 0..5
        std::optional<int>    Presolve;       ///< -1 auto, 0 off, 1 conservative, 2 aggressive
        std::optional<int>    Cuts;           ///< global cut aggressiveness: -1..3
        std::optional<int>    MIPFocus;       ///< 0 default, 1 feasibility, 2 optimality, 3 bound
        std::optional<int>    NumericFocus;   ///< 0..3 numerical-care knob
        std::optional<double> Heuristics;     ///< MIP heuristics effort (0..1)
        std::optional<int>    Threads;        ///< worker thread count; 0 = auto
        std::optional<int>    Seed;           ///< random seed

        // ---- Solution pool (used by mi_opt_multisol) ----
        std::optional<int>    PoolSolutions;
        std::optional<int>    PoolSearchMode;
        std::optional<double> PoolGap;
        std::optional<double> PoolGapAbs;

        // ---- Logging ----
        std::optional<int>         OutputFlag; ///< 0 silent, 1 verbose (Gurobi default 1)
        std::optional<int>         LogToConsole;
        std::optional<std::string> LogFile;

        // ---- Escape hatch for any Gurobi parameter not exposed above ----
        // Keys are the documented Gurobi parameter names, e.g.,
        //   int_params["BarIterLimit"] = 500;
        //   dbl_params["NodefileStart"] = 0.5;
        //   str_params["WLSAccessID"] = "...";
        // See https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html
        std::map<std::string, int>         int_params;
        std::map<std::string, double>      dbl_params;
        std::map<std::string, std::string> str_params;

        // polymorphic copy
        std::unique_ptr<SolverSettings> clone() const override
        {
            return std::make_unique<GurobiSettings>(*this);
        }

        // throws if the Gurobi shared library cannot be dynamically loaded
        void verify_available() const override
        {
            if (!detail::gurobi_available())
            {
                throw std::runtime_error("GurobiSettings: " + detail::gurobi_unavailable_reason());
            }
        }

        /**
         * @brief displays the parameters that have been explicitly set
         */
        std::string print() const
        {
            std::stringstream ss;
            ss << "GurobiSettings (unset fields use Gurobi defaults):\n";
            auto opt_int = [&](const char* n, const std::optional<int>& v)         { if (v) ss << "  " << n << " = " << *v << "\n"; };
            auto opt_dbl = [&](const char* n, const std::optional<double>& v)      { if (v) ss << "  " << n << " = " << *v << "\n"; };
            auto opt_str = [&](const char* n, const std::optional<std::string>& v) { if (v) ss << "  " << n << " = \"" << *v << "\"\n"; };

            opt_dbl("TimeLimit", TimeLimit);
            opt_dbl("WorkLimit", WorkLimit);
            opt_dbl("MemLimit", MemLimit);
            opt_int("SolutionLimit", SolutionLimit);

            opt_dbl("MIPGap", MIPGap);
            opt_dbl("MIPGapAbs", MIPGapAbs);
            opt_dbl("FeasibilityTol", FeasibilityTol);
            opt_dbl("OptimalityTol", OptimalityTol);
            opt_dbl("IntFeasTol", IntFeasTol);

            opt_int("Method", Method);
            opt_int("Presolve", Presolve);
            opt_int("Cuts", Cuts);
            opt_int("MIPFocus", MIPFocus);
            opt_int("NumericFocus", NumericFocus);
            opt_dbl("Heuristics", Heuristics);
            opt_int("Threads", Threads);
            opt_int("Seed", Seed);

            opt_int("PoolSolutions", PoolSolutions);
            opt_int("PoolSearchMode", PoolSearchMode);
            opt_dbl("PoolGap", PoolGap);
            opt_dbl("PoolGapAbs", PoolGapAbs);

            opt_int("OutputFlag", OutputFlag);
            opt_int("LogToConsole", LogToConsole);
            opt_str("LogFile", LogFile);

            for (const auto& kv : int_params) ss << "  [int] "  << kv.first << " = " << kv.second << "\n";
            for (const auto& kv : dbl_params) ss << "  [dbl] "  << kv.first << " = " << kv.second << "\n";
            for (const auto& kv : str_params) ss << "  [str] "  << kv.first << " = \"" << kv.second << "\"\n";
            return ss.str();
        }
    };

    /**
     * @brief Solver-native solution metadata produced by the Gurobi backend.
     *
     * When an OptSolution is produced by a Gurobi solve, OptSolution::external_results
     * points to an instance of this class so the caller can read Gurobi-specific status,
     * node counts, gap, and dual bound. Inspect via dynamic_cast / isinstance.
     */
    struct GurobiSolverResults : ExternalSolverResults
    {
        /// Raw Gurobi status code (see Gurobi's status code table).
        int status = 0;

        /// Simplex / barrier iteration count (Gurobi attribute "IterCount").
        double iter_count = 0.0;

        /// Number of branch-and-bound nodes explored (Gurobi attribute "NodeCount").
        double node_count = 0.0;

        /// Relative MIP optimality gap achieved at termination (Gurobi attribute "MIPGap").
        double mip_gap = 0.0;

        /// Best dual bound found for the MIP (Gurobi attribute "ObjBound").
        double obj_bound = 0.0;

        std::shared_ptr<ExternalSolverResults> clone() const override
        {
            return std::make_shared<GurobiSolverResults>(*this);
        }

        std::string print() const override
        {
            std::stringstream ss;
            ss << "GurobiSolverResults:\n";
            ss << "  status: " << status << "\n";
            ss << "  iter_count: " << iter_count << "\n";
            ss << "  node_count: " << node_count << "\n";
            ss << "  mip_gap: " << mip_gap << "\n";
            ss << "  obj_bound: " << obj_bound << "\n";
            return ss.str();
        }
    };
} // namespace ZonoOpt

#endif
