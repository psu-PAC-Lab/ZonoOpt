#pragma once

/**
 * @file CholeskyUtilities.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Internal utilities for Cholesky factorization using Eigen's LDLT solver.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "Eigen/Sparse"
#include "Eigen/Dense"

#include "Defines.hpp"

namespace ZonoOpt::detail {

    // The completed SimplicialLDLT factorization; stored directly rather than through an extracted
    // copy of L, so factorizing and reading back out don't momentarily double the size of the
    // (possibly much larger) L factor. Not copyable; always held through shared_ptr<const ...>.
    using LDLT_solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>>;

    void affine_set_projection(Eigen::Ref<Eigen::Vector<zono_float, -1>> z, const Eigen::SparseMatrix<zono_float>& A,
        const Eigen::SparseMatrix<zono_float>& AT, const Eigen::Vector<zono_float, -1>& b, const LDLT_solver& AAT_ldlt);

} // end namespace detail
// end namespace ZonoOpt