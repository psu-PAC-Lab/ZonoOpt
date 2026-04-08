#ifndef ZONOOPT_SPARSEMATRIXUTILITIES_HPP_
#define ZONOOPT_SPARSEMATRIXUTILITIES_HPP_

/**
 * @file SparseMatrixUtilities.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Utilities for sparse matrix operations in ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <Eigen/Sparse>
#include <vector>

namespace ZonoOpt::detail
{

    template <typename T>
    Eigen::SparseMatrix<T> hcat(const Eigen::SparseMatrix<T> &A, const Eigen::SparseMatrix<T> &B)
    {
        if (A.rows() != B.rows())
        {
            throw std::invalid_argument("hcat: number of rows must match.");
        }

        Eigen::SparseMatrix<T> C(A.rows(), A.cols() + B.cols());
        std::vector<Eigen::Triplet<T>> tripvec;
        tripvec.reserve(A.nonZeros() + B.nonZeros());

        for (int k=0; k<A.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
            {
                tripvec.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
            }
        }

        for (int k=0; k<B.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(B, k); it; ++it)
            {
                tripvec.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()+A.cols()), it.value());
            }
        }

        C.setFromTriplets(tripvec.begin(), tripvec.end());
        return C;
    }


    template<typename T>
    Eigen::SparseMatrix<T> vcat(const Eigen::SparseMatrix<T> &A, const Eigen::SparseMatrix<T> &B)
    {
        if (A.cols() != B.cols())
        {
            throw std::invalid_argument("vcat: number of columns must match.");
        }

        Eigen::SparseMatrix<T> C(A.rows() + B.rows(), A.cols());
        std::vector<Eigen::Triplet<T>> tripvec;
        tripvec.reserve(A.nonZeros() + B.nonZeros());

        for (int k=0; k<A.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
            {
                tripvec.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
            }
        }

        for (int k=0; k<B.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(B, k); it; ++it)
            {
                tripvec.emplace_back(static_cast<int>(it.row()+A.rows()), static_cast<int>(it.col()), it.value());
            }
        }

        C.setFromTriplets(tripvec.begin(), tripvec.end());
        return C;
    }

    // get triplets for matrix
    template <typename T>
    void get_triplets_offset(const Eigen::SparseMatrix<T> &mat, std::vector<Eigen::Triplet<T>> &triplets,
                const int i_offset, const int j_offset)
    {
        // check validity
        if (i_offset < 0 || j_offset < 0)
        {
            throw std::invalid_argument("get_triplets_offset: offsets must be non-negative.");
        }

        // get triplets
        for (int k=0; k<mat.outerSize(); ++k)
        {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, k); it; ++it)
            {
                triplets.emplace_back(static_cast<int>(it.row() + i_offset), static_cast<int>(it.col() + j_offset), it.value());
            }
        }
    }

    // remove redundant constraints, A*x = b
    template <typename T>
    void remove_redundant_constraints(Eigen::SparseMatrix<T>& A, Eigen::Vector<T,-1>& b)
    {
        // check for empty input matrix
        if (A.rows() == 0 || A.cols() == 0)
        {
            A.resize(0, A.cols());
            b.resize(0);
            return;
        }

        // transpose
        Eigen::SparseMatrix<T> At = A.transpose();

        // compute QR decomposition
        Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> qr;
        qr.analyzePattern(At);
        qr.factorize(At);

        // get the permutation matrix and its indices
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = qr.colsPermutation();
        Eigen::VectorXi P_indices = P.indices();

        // QR solver puts linearly dependent rows at end
        std::vector<Eigen::Triplet<T>> tripvec;
        for (int i=0; i<qr.rank(); i++)
        {
            tripvec.emplace_back(P_indices(i), i, static_cast<T>(1.0));
        }

        Eigen::SparseMatrix<T> P_full (At.cols(), qr.rank());
        P_full.setFromTriplets(tripvec.begin(), tripvec.end());

        // remove redundant constraints
        A = (At * P_full).transpose();
        b = (b.transpose() * P_full).transpose();
    }

    // get triplets for row of row-major sparse matrix
    template <typename T>
    std::vector<Eigen::Triplet<T>> get_triplets_row(const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat, int row)
    {
        std::vector<Eigen::Triplet<T>> triplets;
        if (row < 0 || row >= mat.rows())
        {
            throw std::out_of_range("get_triplets_row: row index out of range.");
        }

        for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(mat, row); it; ++it)
        {
            triplets.emplace_back(row, static_cast<int>(it.col()), it.value());
        }

        return triplets;
    }

    template <typename T>
    struct QREqConsData
    {
        Eigen::SparseMatrix<T> R;
        Eigen::Vector<T, -1> b_tilde;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;
    };

    template <typename T>
    void QREqCons(const Eigen::SparseMatrix<T>& A, const Eigen::Vector<T, -1>& b, QREqConsData<T>& data)
    {
        Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> qr;
        qr.compute(A);
        data.R = qr.matrixR();
        data.b_tilde = qr.matrixQ().transpose()*b;
        data.P = qr.colsPermutation();
    }

} // namespace ZonoOpt::detail
// end namespace ZonoOpt

#endif