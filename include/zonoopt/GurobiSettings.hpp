#ifndef ZONOOPT_GUROBI_SETTINGS_HPP_
#define ZONOOPT_GUROBI_SETTINGS_HPP_

/**
 * @file GurobiSettings.hpp
 * @brief Settings for the dynamically-loaded Gurobi solver backend.
 */

#include <limits>
#include <sstream>
#include <string>

#include "SolverDataStructures.hpp"

namespace ZonoOpt
{
    /**
     * @brief Settings for the Gurobi solver backend.
     *
     * Passing a GurobiSettings instance to a ZonoOpt optimization routine causes the
     * library to attempt dynamic loading of Gurobi and solve via Gurobi's API. If the
     * Gurobi shared library cannot be loaded at runtime, the library silently falls
     * back to the internal ZonoOpt solver with default OptSettings.
     */
    struct GurobiSettings : SolverSettings
    {
        /// display Gurobi output (maps to Gurobi's OutputFlag parameter)
        bool verbose = false;

        /// max wall-clock time for optimization in seconds (Gurobi's TimeLimit parameter)
        double t_max = std::numeric_limits<double>::max();

        /// relative MIP optimality gap (Gurobi's MIPGap parameter)
        double mip_gap = 1e-4;

        /// absolute MIP optimality gap (Gurobi's MIPGapAbs parameter)
        double mip_gap_abs = 1e-10;

        /// number of threads (Gurobi's Threads parameter); 0 = Gurobi's default (typically all cores)
        int threads = 0;

        /// max number of MIP solutions to return in mi_opt_multisol via Gurobi's solution pool
        /// (caps the per-call PoolSolutions; the n_sols argument also caps it)
        int max_pool_solutions = std::numeric_limits<int>::max();

        /// optional path to a Gurobi log file; empty disables file logging
        std::string log_file = "";

        /**
         * @brief displays settings as string
         */
        std::string print() const
        {
            std::stringstream ss;
            ss << "GurobiSettings:" << std::endl;
            ss << "  verbose: " << (verbose ? "true" : "false") << std::endl;
            ss << "  t_max: " << t_max << std::endl;
            ss << "  mip_gap: " << mip_gap << std::endl;
            ss << "  mip_gap_abs: " << mip_gap_abs << std::endl;
            ss << "  threads: " << threads << std::endl;
            ss << "  max_pool_solutions: " << max_pool_solutions << std::endl;
            ss << "  log_file: " << log_file << std::endl;
            return ss.str();
        }
    };
} // namespace ZonoOpt

#endif
