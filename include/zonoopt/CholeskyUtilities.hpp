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

    struct LDLT_data
    {
        Eigen::SparseMatrix<zono_float> L;
        Eigen::DiagonalMatrix<zono_float, -1> Dinv;
        Eigen::PermutationMatrix<-1, -1, int> P, Pinv;
        bool factorized = false;
    };

    void get_LDLT_data(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>>& solver, LDLT_data& data);

    Eigen::Vector<zono_float, -1> solve_LDLT(const LDLT_data& data, const Eigen::Vector<zono_float, -1>& b);

    void affine_set_projection(Eigen::Ref<Eigen::Vector<zono_float, -1>> z, const Eigen::SparseMatrix<zono_float>& A,
        const Eigen::SparseMatrix<zono_float>& AT, const Eigen::Vector<zono_float, -1>& b, const LDLT_data& AAT_ldlt);

} // end namespace detail
// end namespace ZonoOpt