#ifndef ZONOOPT_INTERVALMATRIX_HPP_
#define ZONOOPT_INTERVALMATRIX_HPP_

/**
* @file IntervalMatrix.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief IntervalMatrix class
 * @version 1.0
 * @date 2026-03-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Interval.hpp"
#include "Box.hpp"

namespace ZonoOpt
{
    using namespace detail;

    /**
     * @brief Class for interval matrices (i.e., every element of the matrix is an interval)
     *
     */
    class IntervalMatrix
    {
    public:
        /**
         * @brief IntervalMatrix default constructor
         */
        IntervalMatrix() : rows_(0), cols_(0)
        {
        }

        /**
         * @brief IntervalMatrix constructor using triplets
         * @param rows number of rows
         * @param cols number of columns
         * @param triplets triplets for interval matrix (row, col, interval)
         */
        IntervalMatrix(size_t rows, size_t cols, const std::vector<Eigen::Triplet<Interval>>& triplets);

        /**
         * @brief IntervalMatrix constructor using dense lower and upper bound matrices
         * @param mat_lb lower bound matrix
         * @param mat_ub upper bound matrix
         */
        IntervalMatrix(const Eigen::Matrix<zono_float, -1, -1>& mat_lb,
                       const Eigen::Matrix<zono_float, -1, -1>& mat_ub);

        /**
         * @brief IntervalMatrix constructor from Eigen matrix of intervals
         *
         * @param mat matrix of intervals
         */
        explicit IntervalMatrix(const Eigen::Matrix<Interval, -1, -1>& mat);

        /**
         * @brief Convert to a vector of vectors of Intervals (row-major)
         * 
         * @return vector of vectors of Intervals
         */
        std::vector<std::vector<Interval>> to_array() const;

        /**
         * @brief Convert to triplets
         * 
         * @return std::tuple<int, int, std::vector<Eigen::Triplet<Interval>>> 
         */
        std::tuple<int, int, std::vector<Eigen::Triplet<Interval>>> to_triplets() const;

        /**
         * @brief Get center matrix
         * @return center matrix
         *
         * Each element of center matrix is the center of the corresponding interval in the interval matrix
         */
        Eigen::SparseMatrix<zono_float> center() const;

        /**
         * @brief Get diameter matrix
         * @return diameter matrix
         *
         * Each element of the diameter matrix is the width of the corresponding interval in the interval matrix
         */
        Eigen::SparseMatrix<zono_float> diam() const;

        /**
         * @brief Get radius matrix
         * @return radius matrix
         *
         * Returns the IntervalMatrix with each interval shifted to be centered at zero
         */
        IntervalMatrix radius() const;

        /**
         * @brief Get width of interval matrix
         * @return width of interval matrix
         *
         * Specifically, this returns the max width for any interval in the interval matrix
         */
        zono_float width() const;

        /**
         * @brief Interval matrix intersection
         * 
         * @param other interval matrix to intersect with
         * @return IntervalMatrix 
         */
        IntervalMatrix intersect(const IntervalMatrix& other) const;

        /**
         * @brief Interval matrix interval hull
         * 
         * @param other interval matrix to compute hull with
         * @return IntervalMatrix 
         */
        IntervalMatrix interval_hull(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix multiplication with vector
         * @param v rhs vector
         * @return resulting box
         */
        Box operator*(const Eigen::Vector<zono_float, -1>& v) const;

        /**
         * @brief IntervalMatrix multiplication with Box
         * @param box rhs box
         * @return resulting box
         */
        Box operator*(const Box& box) const;

        /**
         * @brief IntervalMatrix scalar multiplication
         * @param alpha scalar multiplier
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(zono_float alpha) const;

        /**
         * @brief IntervalMatrix scalar multiplication in-place
         * @param alpha scalar multiplier
         */
        void operator*=(zono_float alpha);

        /**
         * @brief IntervalMatrix scalar multiplication
         * @param alpha scalar multiplier
         * @param A interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator*(zono_float alpha, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix elementwise multiplication by interval
         * @param interval interval to multiply
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const Interval& interval) const;

        /**
         * @brief IntervalMatrix elementwise multiplication by interval
         * @param interval interval to multiply
         * @param A interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator*(const Interval& interval, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix elementwise multiplication by interval in-place
         * @param interval interval to multiply
         */
        void operator*=(const Interval& interval);

        /**
         * @brief IntervalMatrix elementwise division by scalar
         * @param alpha scalar to divide
         * @return resulting interval matrix (this / alpha)
         */
        IntervalMatrix operator/(zono_float alpha) const;

        /**
         * @brief IntervalMatrix elementwise scalar division
         * @param alpha scalar
         * @param A matrix to divide
         * @return resulting interval matrix (alpha / this)
         */
        friend IntervalMatrix operator/(zono_float alpha, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix elementwise division by scalar
         * @param alpha scalar to divide
         */
        void operator/=(zono_float alpha);

        /**
         * @brief IntervalMatrix elementwise division by interval
         * @param interval interval to divide
         * @return resulting interval matrix (this / interval)
         */
        IntervalMatrix operator/(const Interval& interval) const;

        /**
         * @brief IntervalMatrix elementwise interval division
         * @param interval interval
         * @param A matrix to divide
         * @return resulting interval matrix (interval / this)
         */
        friend IntervalMatrix operator/(const Interval& interval, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix elementwise division by interval in-place
         * @param interval interval to divide
         */
        void operator/=(const Interval& interval);

        /**
         * @brief IntervalMatrix multiplication with sparse matrix
         * @param A rhs matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const;

        /**
         * @brief IntervalMatrix multiplication with sparse matrix
         * @param A lhs matrix
         * @param B rhs interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const IntervalMatrix& B);

        /**
         * @brief IntervalMatrix multiplication with dense matrix
         * @param A rhs matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const Eigen::Matrix<zono_float, -1, -1>& A) const;

        /**
         * @brief IntervalMatrix multiplication with dense matrix
         * @param A lhs matrix
         * @param B rhs interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator*(const Eigen::Matrix<zono_float, -1, -1>& A, const IntervalMatrix& B);

        /**
         * @brief IntervalMatrix multiplication with another IntervalMatrix
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix multiplication with another IntervalMatrix in-place
         * @param other rhs interval matrix
         */
        void operator*=(const IntervalMatrix& other);

        /**
         * @brief IntervalMatrix addition
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator+(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix addition in-place
         * @param other rhs interval matrix
         */
        void operator+=(const IntervalMatrix& other);

        /**
         * @brief IntervalMatrix element-wise addition with interval
         * @param interval interval to add
         * @return resulting interval matrix
         */
        IntervalMatrix operator+(const Interval& interval) const;

        /**
         * @brief IntervalMatrix element-wise addition with interval
         * @param interval interval to add
         * @param mat interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator+(const Interval& interval, const IntervalMatrix& mat);

        /**
         * @brief IntervalMatrix element-wise addition with interval in-place
         * @param interval interval to add
         */
        void operator+=(const Interval& interval);

        /**
         * @brief IntervalMatrix element-wise addition with scalar
         * @param alpha scalar to add
         * @return resulting interval matrix
         */
        IntervalMatrix operator+(zono_float alpha) const;

        /**
         * @brief IntervalMatrix element-wise addition with scalar
         * @param alpha scalar to add
         * @param A interval matrix
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator+(zono_float alpha, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix element-wise addition with scalar in-place
         * @param alpha scalar to add
         */
        void operator+=(zono_float alpha);


        /**
         * @brief IntervalMatrix subtraction
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator-(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix subtraction in-place
         * @param other rhs interval matrix
         */
        void operator-=(const IntervalMatrix& other);


        /**
         * @brief IntervalMatrix elementwise interval subtraction
         * @param interval interval to subtract
         * @return resulting interval matrix
         */
        IntervalMatrix operator-(const Interval& interval) const;

        /**
         * @brief IntervalMatrix elementwise interval subtraction
         * @param interval interval
         * @param mat interval matrix to subtract
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator-(const Interval& interval, const IntervalMatrix& mat);

        /**
         * @brief IntervalMatrix elementwise interval subtraction in-place
         * @param interval interval
         */
        void operator-=(const Interval& interval);

        /**
         * @brief IntervalMatrix elementwise scalar subtraction
         * @param alpha scalar to subtract
         * @return resulting interval matrix
         */
        IntervalMatrix operator-(zono_float alpha) const;

        /**
         * @brief IntervalMatrix elementwise scalar subtraction
         * @param alpha scalar
         * @param A interval matrix to subtract
         * @return resulting interval matrix
         */
        friend IntervalMatrix operator-(zono_float alpha, const IntervalMatrix& A);

        /**
         * @brief IntervalMatrix elementwise scalar subtraction in-place
         * @param alpha scalar to subtract
         */
        void operator-=(zono_float alpha);

        /**
         * @brief IntervalMatrix negation
         * @return negated interval matrix
         */
        IntervalMatrix operator-() const;

        /**
         * @brief Interval matrix intersection oeprator
         * 
         * @param other interval matrix to intersect with
         * @return IntervalMatrix 
         */
        IntervalMatrix operator&(const IntervalMatrix& other) const;

        /**
         * @brief Interval matrix interval hull operator
         * 
         * @param other interval matrix to compute hull with
         * @return IntervalMatrix 
         */
        IntervalMatrix operator|(const IntervalMatrix& other) const;

        /**
         * @brief Get number of rows
         * @return number of rows
         */
        size_t rows() const { return rows_; }

        /**
         * @brief Get number of columns
         * @return number of cols
         */
        size_t cols() const { return cols_; }

        /**
         * @brief Print method
         * @return string display of IntervalMatrix
         */
        std::string print() const;

        /**
         * @brief print to ostream
         * @param os ostream
         * @param interval_matrix reference to interval matrix
         * @return ostream
         */
        friend std::ostream& operator<<(std::ostream& os, const IntervalMatrix& interval_matrix);

    private:
        size_t rows_, cols_;
        std::vector<std::vector<std::pair<size_t, Interval>>> mat_; // rows->cols->vals
    };

}


#endif