#ifndef ZONOOPT_SCIP_SOLVER_HPP_
#define ZONOOPT_SCIP_SOLVER_HPP_

/**
 * @file SCIPSolver.hpp
 * @brief C++ interface to the SCIP optimization solver via dynamic loading.
 * 
 */

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "SolverDataStructures.hpp"
#include "SCIPSettings.hpp"

namespace ZonoOpt
{
namespace detail
{

/**
 * @brief True if the SCIP shared library could be dynamically loaded at program startup.
 */
bool scip_available();

/**
 * @brief Solve the QP that backs ConZono::qp_opt using SCIP.
 *
 *   min  0.5 * xi^T P xi + q^T xi + c
 *   s.t. A xi = b,  xi_lb <= xi <= xi_ub
 *
 * If P is the zero matrix the problem is solved as an LP. Otherwise the quadratic
 * objective is modeled via the standard epigraph reformulation: add auxiliary
 * variable t and the quadratic constraint 0.5 * xi^T P xi - t <= 0, then minimize
 * q^T xi + t (constant c is added back to OptSolution::J).
 *
 * @throws std::runtime_error if the SCIP API is unavailable or a fatal SCIP API error occurs during model construction or solve.
 */
OptSolution solve_qp_scip(const Eigen::SparseMatrix<zono_float>& P,
                          const Eigen::Vector<zono_float, -1>& q,
                          zono_float c,
                          const Eigen::SparseMatrix<zono_float>& A,
                          const Eigen::Vector<zono_float, -1>& b,
                          const Eigen::Vector<zono_float, -1>& xi_lb,
                          const Eigen::Vector<zono_float, -1>& xi_ub,
                          const SCIPSettings& settings);

/**
 * @brief Solve the MIQP that backs HybZono::mi_opt using SCIP.
 *
 * Identical formulation to solve_qp_scip but variables in [bin_start, bin_start+bin_count)
 * are restricted to {0, 1}. When zero_one_form is false, the substitution
 * xi_i = 2*y_i - 1 is applied internally and undone after solving.
 *
 * @throws std::runtime_error if the SCIP API is unavailable or a fatal SCIP API error occurs during model construction or solve.
 */
OptSolution solve_miqp_scip(const Eigen::SparseMatrix<zono_float>& P,
                            const Eigen::Vector<zono_float, -1>& q,
                            zono_float c,
                            const Eigen::SparseMatrix<zono_float>& A,
                            const Eigen::Vector<zono_float, -1>& b,
                            const Eigen::Vector<zono_float, -1>& xi_lb,
                            const Eigen::Vector<zono_float, -1>& xi_ub,
                            int bin_start, int bin_count,
                            bool zero_one_form,
                            const SCIPSettings& settings);

/**
 * @brief Return up to n_sols feasible MIQP solutions collected by SCIP's solution storage.
 *
 * SCIP collects every feasible solution found during branch-and-bound into its solution
 * storage (capped by "limits/maxsol", which we set to n_sols). After solving, we iterate
 * the storage and return one OptSolution per stored solution, ordered best-first by
 * SCIP's internal ranking. Note: these are "feasible solutions found during search",
 * not a "find n best" enumeration — same semantics as ZonoOpt's internal multisol.
 *
 * @throws std::runtime_error if the SCIP API is unavailable or a fatal SCIP API error occurs during model construction or solve.
 */
std::vector<OptSolution> solve_miqp_scip_multisol(const Eigen::SparseMatrix<zono_float>& P,
                                                  const Eigen::Vector<zono_float, -1>& q,
                                                  zono_float c,
                                                  const Eigen::SparseMatrix<zono_float>& A,
                                                  const Eigen::Vector<zono_float, -1>& b,
                                                  const Eigen::Vector<zono_float, -1>& xi_lb,
                                                  const Eigen::Vector<zono_float, -1>& xi_ub,
                                                  int bin_start, int bin_count,
                                                  bool zero_one_form,
                                                  int n_sols,
                                                  const SCIPSettings& settings);

} // namespace detail
} // namespace ZonoOpt

#endif
