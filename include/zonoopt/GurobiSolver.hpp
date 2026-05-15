#ifndef ZONOOPT_GUROBI_SOLVER_HPP_
#define ZONOOPT_GUROBI_SOLVER_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "SolverDataStructures.hpp"

namespace ZonoOpt
{
namespace detail
{

/**
 * @brief Returns true if Gurobi was successfully loaded at program startup.
 *
 * Safe to call regardless of build configuration; if Gurobi is not installed
 * or could not be loaded, returns false and the caller should fall back to
 * the internal solver.
 */
bool gurobi_available();

/**
 * @brief Solve the QP that backs ConZono::qp_opt using Gurobi.
 *
 * Problem solved:
 *   min  0.5 * xi^T P xi + q^T xi + c
 *   s.t. A xi = b
 *        xi_lb <= xi <= xi_ub
 *
 * Returns an OptSolution with z, J, run_time, converged, and infeasible populated.
 * Throws std::runtime_error if Gurobi is loaded but encounters a fatal API error.
 */
OptSolution solve_qp_gurobi(const Eigen::SparseMatrix<zono_float>& P,
                            const Eigen::Vector<zono_float, -1>& q,
                            zono_float c,
                            const Eigen::SparseMatrix<zono_float>& A,
                            const Eigen::Vector<zono_float, -1>& b,
                            const Eigen::Vector<zono_float, -1>& xi_lb,
                            const Eigen::Vector<zono_float, -1>& xi_ub,
                            const OptSettings& settings);

/**
 * @brief Solve the MIQP that backs HybZono::mi_opt using Gurobi.
 *
 * Problem solved (in the variable xi):
 *   min  0.5 * xi^T P xi + q^T xi + c
 *   s.t. A xi = b
 *        xi_lb <= xi <= xi_ub
 *        xi_i in {xi_lb_i, xi_ub_i} for i in [bin_start, bin_start + bin_count)
 *
 * If zero_one_form is true, binary variables are in {0, 1}.
 * If zero_one_form is false, binary variables are in {-1, 1} and the solver
 * internally substitutes y_i = (xi_i + 1) / 2 so Gurobi can use type='B' (which
 * requires {0,1}); the recovered xi is returned via OptSolution.z.
 *
 * Returns an OptSolution with z, J, run_time, converged, and infeasible populated.
 */
OptSolution solve_miqp_gurobi(const Eigen::SparseMatrix<zono_float>& P,
                              const Eigen::Vector<zono_float, -1>& q,
                              zono_float c,
                              const Eigen::SparseMatrix<zono_float>& A,
                              const Eigen::Vector<zono_float, -1>& b,
                              const Eigen::Vector<zono_float, -1>& xi_lb,
                              const Eigen::Vector<zono_float, -1>& xi_ub,
                              int bin_start, int bin_count,
                              bool zero_one_form,
                              const OptSettings& settings);

/**
 * @brief Enumerate up to n_sols MIQP solutions using Gurobi's solution pool.
 *
 * Uses PoolSearchMode=2 (find n best solutions) with PoolSolutions=n_sols.
 * Returns a vector of OptSolution, one per pool entry (size up to n_sols, possibly 0
 * if the problem is infeasible).
 */
std::vector<OptSolution> solve_miqp_gurobi_multisol(const Eigen::SparseMatrix<zono_float>& P,
                                                    const Eigen::Vector<zono_float, -1>& q,
                                                    zono_float c,
                                                    const Eigen::SparseMatrix<zono_float>& A,
                                                    const Eigen::Vector<zono_float, -1>& b,
                                                    const Eigen::Vector<zono_float, -1>& xi_lb,
                                                    const Eigen::Vector<zono_float, -1>& xi_ub,
                                                    int bin_start, int bin_count,
                                                    bool zero_one_form,
                                                    int n_sols,
                                                    const OptSettings& settings);

} // namespace detail
} // namespace ZonoOpt

#endif
