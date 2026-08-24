#pragma once

/**
 * @file Box.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Box classe
 * @version 1.0
 * @date 2026-03-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <vector>
#include <set>
#include <utility>
#include <ostream>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Interval.hpp"
#include "Defines.hpp"

namespace ZonoOpt
{
    using namespace detail;

    namespace detail
    {
        class MI_Box;
    }

    /**
     * @brief Box (i.e., interval vector) class
     *
     * @warning Unlike HybZono and its subclasses, Box operators have interval-arithmetic
     * semantics, not set-operation semantics. Specifically, Box::operator* is elementwise
     * interval multiplication (not Cartesian product) and Box::operator- is interval
     * subtraction [a,b] - [c,d] = [a-d, b-c] (not the Pontryagin difference [a-c, b-d]).
     * Use the named methods cartesian_product() and pontry_diff() for the set-operation
     * meanings. Box::operator+ is elementwise interval addition, which for boxes coincides
     * exactly with the Minkowski sum; minkowski_sum() is a named alias for it.
     */
    class Box
    {
    public:
        // constructors

        /**
         * @brief Default constructor
         */
        Box() = default;

        /**
         * @brief Default construct with size specified
         * @param size dimension of box
         */
        explicit Box(const size_t size);

        /**
         * @brief Constructor using vector of intervals
         * @param vals vector of intervals
         */
        explicit Box(const std::vector<Interval>& vals);

        /**
         * @brief Constructor using Eigen vector of intervals
         * @param vals Eigen vector of intervals
         */
        explicit Box(const Eigen::Vector<Interval, -1>& vals);

        /**
         * @brief Constructor from intervals of lower and upper bounds
         * @param x_lb vector of lower bounds
         * @param x_ub vector of upper bounds
         * @throws std::invalid_argument if x_lb and x_ub do not have the same size.
         */
        Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub);

        /**
         * @brief Convert to vector of intervals
         * 
         * @return std::vector<Interval> 
         */
        std::vector<Interval> to_array() const;

        // virtual destructor
        /**
         * @brief Virtual destructor
         */
        virtual ~Box() = default;

        // copy assignment
        /**
         * @brief Copy assignment
         * @param other other Box object
         * @return this = other
         */
        Box& operator=(const Box& other);

        // copy constructor
        /**
         * @brief Copy constructor
         * @param other other Box object
         */
        Box(const Box& other);

        // element-wise assignment, access

        /**
         * @brief Element-wise access
         * @param i index
         * @return Interval for element i in Box
         * @throws std::out_of_range if i is outside [0, size()).
         */
        Interval get_element(int i) const;


        /**
         * @brief Element-wise assignment
         * @param i index
         * @param val interval to assign
         */
        void set_element(int i, const Interval& val);

        /**
         * @brief Get size of Box object
         * @return size of box
         */
        size_t size() const;

        // project onto box

        /**
         * @brief Projects vector onto the Box
         * @param x vector reference
         * @throws std::invalid_argument if x does not have the same size as the Box.
         */
        virtual void project(Eigen::Ref<Eigen::Vector<zono_float, -1>> x) const;

        /**
         * @brief Clone operation
         * @return pointer to newly created object
         */
        virtual Box* clone() const;

        // access bounds
        /**
         * @brief Get lower bounds
         * @return const reference to lower bounds
         */
        const Eigen::Vector<zono_float, -1>& lower() const { return x_lb; }

        /**
         * @brief Get upper bounds
         * @return const reference to upper bounds
         */
        const Eigen::Vector<zono_float, -1>& upper() const { return x_ub; }

        // width of box
        /**
         * @brief Get width of box
         * @return width of box
         *
         * Specifically, this returns the max width for any interval in the box
         */
        zono_float width() const;

        /**
         * @brief Get radius of box
         * @return radius of box
         *
         * Returns box with intervals centered at zero with width equal to the width of the original box
         */
        Box radius() const;

        /**
         * @brief Get center of box
         * @return center of box
         */
        Eigen::Vector<zono_float, -1> center() const;

        /**
         * @brief Determine whether box is empty (any contained interval is empty)
         * @return true if box is empty
         */
        bool is_empty() const;

        /**
         * @brief Determine whether box is single-valued (i.e., all intervals have width 0 within numerical tolerance)
         * @return true if box is single-valued
         */
        bool is_single_valued() const;

        /**
         * @brief Box intersection
         *
         * @param other other
         * @return intersection of this and other
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box intersect(const Box& other) const;

        /**
         * @brief Box interval hull
         *
         * @param other other box
         * @return interval hull of this and other
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box interval_hull(const Box& other) const;

        /**
         * @brief Box intersection over a subset of dimensions
         *
         * @param other other box, must have size equal to dims.size()
         * @param dims dimensions of this box over which the intersection is computed
         * @return intersection of this and other over the specified dimensions
         * @throws std::invalid_argument if dims.size() does not equal other.size(), or if any entry in dims is not a valid dimension of this box.
         *
         * This computes the exact intersection {x in this : x(dims) in other}. Repeated
         * indices in dims are intersected cumulatively.
         */
        Box intersection_over_dims(const Box& other, const std::vector<int>& dims) const;

        /**
         * @brief Check vector continment
         *
         * @param v vector
         * @return true if v is contained in box, false otherwise
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        bool contains(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Check set containment
         *
         * @param other box to check if subset
         * @return true if other is a subset of this
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        bool contains_set(const Box& other) const;

        /**
         * @brief Cartesian product of this box and other
         *
         * @param other other box
         * @return Box of dimension size() + other.size() with this box's intervals first
         *
         * @note Box::operator* is elementwise interval multiplication, NOT the Cartesian
         * product. This differs from HybZono, where operator* is the Cartesian product.
         */
        Box cartesian_product(const Box& other) const;

        /**
         * @brief Minkowski sum of this box and other
         *
         * @param other other box
         * @return Minkowski sum of this and other
         *
         * @note For boxes, the Minkowski sum coincides exactly with elementwise interval
         * addition, so this is a named alias for Box::operator+.
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box minkowski_sum(const Box& other) const;

        /**
         * @brief Pontryagin difference of this box and other, i.e., {x : x + y is in this for all y in other}
         *
         * @param other subtrahend
         * @return Pontryagin difference of this and other
         *
         * @warning Box::operator- is interval subtraction, i.e., [a,b] - [c,d] = [a-d, b-c],
         * which is NOT the Pontryagin difference [a-c, b-d]. This differs from HybZono,
         * where operator- is the Pontryagin difference.
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box pontry_diff(const Box& other) const;

        /**
         * @brief Projects the box onto the dimensions specified in dims
         *
         * @param dims vector of dimensions
         * @return Box of dimension dims.size() whose k-th interval is this box's dims[k]-th interval
         *
         * dims need not be sorted and may contain repeats.
         *
         * @throws std::invalid_argument if any entry in dims is not a valid dimension of the Box.
         */
        Box project_onto_dims(const std::vector<int>& dims) const;

        // operator overloading

        /**
         * @brief Elementwise addition
         * @param other rhs box
         * @return enclosure of this + other (elementwise)
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator+(const Box& other) const;

        /**
         * @brief Elementwise addition in-place
         * @param other other box
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        void operator+=(const Box& other);

        /**
         * @brief Elementwise addition with vector
         * @param v vector to add
         * @param box box to add
         * @return enclosure of v + box (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        Box operator+(const Eigen::Vector<zono_float, -1>& v) const;


        /**
         * @brief Elementwise addition with vector in-place
         * @param v vector
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        void operator+=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise addition with vector
         * @param v vector to add
         * @param box box to add
         * @return enclosure of v + box (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        friend Box operator+(const Eigen::Vector<zono_float, -1>& v, const Box& box);

        /**
         * @brief Elementwise subtraction
         * @param other rhs box
         * @return enclosure of this - other (elementwise)
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator-(const Box& other) const;

        /**
         * @brief Elementwise subtraction in-place
         * @param other other box
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        void operator-=(const Box& other);

        /**
         * @brief Elementwise subtraction with vector
         * @param v vector to subtract
         * @return enclosure of this - v (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        Box operator-(const Eigen::Vector<zono_float, -1>& v) const;

        /**
         * @brief Elementwise subtraction with vector in-place
         * @param v vector to subtract
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        void operator-=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise subtraction with vector
         * @param v vector
         * @param box box to subtract
         * @return enclosure of v - box (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        friend Box operator-(const Eigen::Vector<zono_float, -1>& v, const Box& box);

        /**
         * @brief Elementwise multiplication
         * @param other rhs box
         * @return enclosure of this * other (elementwise)
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator*(const Box& other) const;

        /**
         * @brief Elementwise multiplication in-place
         * @param other other box
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        void operator*=(const Box& other);

        /**
         * @brief Elementwise multiplication with scalar
         * @param alpha scalar multiplier
         * @return enclosure of alpha * this (elementwise)
         */
        Box operator*(zono_float alpha) const;

        /**
         * @brief Elementwise multiplication with scalar in-place
         * @param alpha scalar multiplier
         */
        void operator*=(zono_float alpha);

        /**
         * @brief Elementwise multiplication with scalar
         * @param alpha scalar multiplier
         * @param box box to multiply
         * @return enclosure of alpha * box (elementwise)
         */
        friend Box operator*(zono_float alpha, const Box& box);

        /**
         * @brief Elementwise multiplication with vector
         * @param v vector to multiply
         * @return enclosure of this * v (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        Box operator*(const Eigen::Vector<zono_float, -1>& v) const;

        /**
         * @brief Elementwise multiplication with vector in-place
         * @param v vector to multiply
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        void operator*=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise multiplication with vector
         * @param v vector to multiply
         * @param box box to multiply
         * @return enclosure of v * box (elementwise)
         * @throws std::invalid_argument if v does not have the same size as the Box.
         */
        friend Box operator*(const Eigen::Vector<zono_float, -1>& v, const Box& box);

        /**
         * @brief Elementwise multiplication with interval
         * @param interval interval to multiply
         * @return enclosure of this * interval (elementwise)
         */
        Box operator*(const Interval& interval) const;

        /**
         * @brief Elementwise multiplication with interval in-place
         * @param interval interval to multiply
         */
        void operator*=(const Interval& interval);

        /**
         * @brief Elementwise multiplication with interval
         * @param interval interval to multiply
         * @param box box to multiply
         * @return enclosure of interval * box (elementwise)
         */
        friend Box operator*(const Interval& interval, const Box& box);

        /**
         * @brief Linear map
         * @param A matrix to multiply
         * @return enclosure of A * this
         */
        friend Box operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const Box& box);

        /**
         * @brief Linear map
         * @param A matrix to multiply
         * @return enclosure of A * this
         */
        friend Box operator*(const Eigen::Matrix<zono_float, -1, -1>& A, const Box& box);

        /**
         * @brief Elementwise division
         * @param other rhs box
         * @return enclosure of this / other (elementwise)
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator/(const Box& other) const;

        /**
         * @brief Elementwise division in-place
         * @param other rhs box
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        void operator/=(const Box& other);

        /**
         * @brief Elementwise division with scalar
         * @param alpha scalar divisor
         * @return enclosure of this / alpha (elementwise)
         */
        Box operator/(zono_float alpha) const;

        /**
         * @brief Elementwise division with scalar
         * @param alpha scalar divisor
         */
        void operator/=(zono_float alpha);

        /**
         * @brief Elementwise division with scalar
         * @param alpha scalar divisor
         * @param box box to divide
         * @return enclosure of alpha / box (elementwise)
         */
        friend Box operator/(zono_float alpha, const Box& box);

        /**
         * @brief Elementwise division with interval
         * @param interval interval to divide
         * @return enclosure of this / interval (elementwise)
         */
        Box operator/(const Interval& interval) const;

        /**
         * @brief Elementwise division with interval in-place
         * @param interval interval to divide
         */
        void operator/=(const Interval& interval);

        /**
         * @brief Elementwise division with interval
         * @param interval interval
         * @param box box to divide
         * @return enclosure of interval / box (elementwise)
         */
        friend Box operator/(const Interval& interval, const Box& box);

        /**
         * @brief Unary minus: returns -1 * this
         * @return enclosure of -this
         */
        Box operator-() const;

        /**
         * @brief Box intersection
         *
         * @param other other box
         * @return Intersection of this and other
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator&(const Box& other) const;

        /**
         * @brief Box interval hull
         *
         * @param other box
         * @return Interval hull of this and other
         * @throws std::invalid_argument if this and other have inconsistent dimensions.
         */
        Box operator|(const Box& other) const;

        /**
         * @brief Set containment operator
         * 
         * @param other 
         * @return true if this is subset of other
         */
        bool operator<=(const Box& other) const;

        /**
         * @brief Set containment operator
         * 
         * @param other 
         * @return true if this is superset of other
         */
        bool operator>=(const Box& other) const;

        /**
         * @brief Set equality operator
         * 
         * @param other 
         * @return true if sets are equal
         */
        bool operator==(const Box& other) const;

        // interval contractors

        /**
         * @brief Interval contractor
         * @param A constraint matrix (row-major)
         * @param b constraint vector
         * @param iter number of contractor iterations
         * @return flag indicating that the contractor did not detect that A*x=b and the box do not intersect
         *
         * Executes a forward-backward interval contractor for the equality constraint A*x=b.
         * For points x in the box, this shrinks the box without removing any points x that satisfy A*x=b.
         * If the contractor detects that the box does not intersect A*x=b, then this function will return false.
         *
         * @throws std::invalid_argument if iter is not positive.
         */
        bool contract(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const Eigen::Vector<zono_float, -1>& b,
                      int iter);

        /**
         * @brief Interval contractor over a subset of the dimensions of the box
         * @param A_rm constraint matrix, row major
         * @param b constraint vector
         * @param iter number of contractor iterations
         * @param A constraint matrix, column major
         * @param inds box dimension indices
         * @param tree_search_depth how deep to search constraint tree
         * @return flag indicating that the contractor did not detect that A*x=b and the box do not intersect
         *
         * This is a forward-backward contractor over a subset of the dimensions of the box.
         * This detects what other dimensions are affected up to a specified search depth prior to executing the contractor.
         *
         * @throws std::invalid_argument if iter is not positive or if A_rm does not match A.
         */
        bool contract_subset(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                             const Eigen::Vector<zono_float, -1>& b, int iter,
                             const Eigen::SparseMatrix<zono_float>& A, const std::set<int>& inds,
                             int tree_search_depth);

        /**
         * @brief Linear map of box based on interval arithmetic
         * @param A map matrix (dense)
         * @return Linear mapped box
         * @throws std::invalid_argument if A's number of columns does not equal the size of the Box.
         */
        Box linear_map(const Eigen::Matrix<zono_float, -1, -1>& A) const;

        /**
         * @brief Linear map of box based on interval arithmetic
         * @param A map matrix (sparse row major)
         * @return Linear mapped box
         * @throws std::invalid_argument if A's number of columns does not equal the size of the Box.
         */
        Box linear_map(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const;

        /**
         * @brief Linear map of box based on interval arithmetic
         * @param A map matrix (sparse column major)
         * @return Linear mapped box
         * @throws std::invalid_argument if A's number of columns does not equal the size of the Box.
         */
        Box linear_map(const Eigen::SparseMatrix<zono_float>& A) const;

        /**
         * @brief Affine map R*[x] + s of the box based on interval arithmetic
         * @param R map matrix (dense)
         * @param s vector offset, defaults to zero
         * @return affine mapped box
         * @throws std::invalid_argument if R, s, and the Box have inconsistent dimensions.
         */
        Box affine_map(const Eigen::Matrix<zono_float, -1, -1>& R,
                       const Eigen::Vector<zono_float, -1>& s = Eigen::Vector<zono_float, -1>()) const;

        /**
         * @brief Affine map R*[x] + s of the box based on interval arithmetic
         * @param R map matrix (sparse row major)
         * @param s vector offset, defaults to zero
         * @return affine mapped box
         * @throws std::invalid_argument if R, s, and the Box have inconsistent dimensions.
         */
        Box affine_map(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& R,
                       const Eigen::Vector<zono_float, -1>& s = Eigen::Vector<zono_float, -1>()) const;

        /**
         * @brief Affine map R*[x] + s of the box based on interval arithmetic
         * @param R map matrix (sparse column major)
         * @param s vector offset, defaults to zero
         * @return affine mapped box
         * @throws std::invalid_argument if R, s, and the Box have inconsistent dimensions.
         */
        Box affine_map(const Eigen::SparseMatrix<zono_float>& R,
                       const Eigen::Vector<zono_float, -1>& s = Eigen::Vector<zono_float, -1>()) const;

        /**
         * @brief Linear map with vector
         * @param x vector
         * @return Interval
         * @throws std::invalid_argument if x does not have the same size as the Box.
         */
        Interval dot(const Eigen::Vector<zono_float, -1>& x) const;

        /**
         * @brief Computes the support function of the box in the direction d
         *
         * @param d vector defining direction for support function
         * @return support
         *
         * Computes max_{x in [x]} <x, d>. Unlike HybZono::support this is exact and
         * closed-form (no solver settings are required), since each variable appears
         * exactly once in the sum and there is no dependency problem.
         *
         * @throws std::invalid_argument if d does not have the same size as the Box.
         */
        zono_float support(const Eigen::Vector<zono_float, -1>& d) const;

        /**
         * @brief Permutes in place using permutation matrix, i.e., [x] <- P*[x]
         * @param P permutation matrix
         */
        void permute(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P);

        /**
         * @brief Print method
         * @return string display of Box
         */
        std::string print() const;

        /**
         * @brief print to ostream
         * @param os ostream
         * @param box reference to box
         * @return ostream
         */
        friend std::ostream& operator<<(std::ostream& os, const Box& box);

    protected:
        // members
        /// vector of lower bounds
        Eigen::Vector<zono_float, -1> x_lb;

        /// vector of upper bounds
        Eigen::Vector<zono_float, -1> x_ub;

    private:
        // back end for contraction operator

        /**
         * @brief Back-end for contractor
         * @param A constraint matrix (row-major)
         * @param b constraint vector
         * @param iter number of contractor iterations
         * @param constraints constraints to consider
         * @return flag indicating that the contractor did not detect that A*x=b and the box do not intersect
         */
        virtual bool contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                                     const Eigen::Vector<zono_float, -1>& b, const int iter,
                                     const std::set<int>& constraints);

        /**
         * @brief Forward propagation step of contractor
         * @param A constraint matrix (row-major)
         * @param b constraint vector
         * @param constraints constraints to consider
         * @return success flag indicating that the contractor did not detect that A*x=b and the box do not intersect
         */
        bool contract_forward(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                              const Eigen::Vector<zono_float, -1>& b, const std::set<int>& constraints);

        /**
         * @brief Backward propagation step of contractor
         * @param A constraint matrix (row-major)
         * @param b constraint vector
         * @param constraints constraints to consider
         */
        void contract_backward(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                               const Eigen::Vector<zono_float, -1>& b, const std::set<int>& constraints);

        /**
         * @brief Search constraint tree for affected dimensions of box (recursive)
         * @param A constraint matrix (column major)
         * @param A_rm constraint matrix (row major)
         * @param constraints constraints to search
         * @param vars indices to consider
         * @param new_vars all affected indices (reference)
         * @param depth current depth in constraint tree
         * @param max_depth max depth to search constraint tree
         */
        static void get_vars_cons(const Eigen::SparseMatrix<zono_float>& A,
                                  const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                                  std::set<int>& constraints, std::set<int>& vars, const std::set<int>& new_vars,
                                  int depth, int max_depth);

    private:
        friend class detail::MI_Box;
    };

    namespace detail 
    {

        // mixed-integer box (some bounds are fixed)
        /**
        * @brief Mixed-integer box
        *
        * Extends Box class to include variables for which may only take their upper or lower bound
        * (they may not take any value on the interior).
        *
        * @note The set-operation methods inherited from Box (cartesian_product, minkowski_sum,
        * pontry_diff, project_onto_dims, affine_map) are not virtual and always return a plain
        * Box computed from the interval relaxation; integrality information is discarded.
        */
        class MI_Box final : public Box
        {
        public:
            // constructors

            /**
            * @brief default constructor
            */
            MI_Box() = default;

            /**
            * @brief Constructor for MI_Box
            * @param x_lb vector of lower bounds
            * @param x_ub vector of upper bounds
            * @param idx_b indices of binary variables {start index, number of binaries}
            * @param zero_one_form flag indicating whether binary variables are in {0,1} form (true) or {-1,1} form (false)
            * @throws std::invalid_argument if x_lb and x_ub do not have the same size.
            * @throws std::out_of_range if the binary indices in idx_b are outside [0, size()).
            */
            MI_Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub,
                const std::pair<int, int>& idx_b, bool zero_one_form);

            /**
            * @brief Constructor for MI_Box from vector of intervals
            * @param intervals vector intervals
            * @param idx_b indices of binary variables {start index, number of binaries}
            * @param zero_one_form flag indicating whether binary variables are in {0,1} form (true) or {-1,1} form (false)
            * @throws std::out_of_range if the binary indices in idx_b are outside [0, size()).
            */
            MI_Box(const std::vector<Interval>& intervals, const std::pair<int, int>& idx_b, bool zero_one_form);

            // clone operation
            Box* clone() const override;

            // project
            void project(Eigen::Ref<Eigen::Vector<zono_float, -1>> x) const override;

            // get binary indices
            /**
            * @brief Get binary indices
            * @return reference to binary indices {start index, number of binaries}
            */
            const std::pair<int, int>& binary_indices() const { return idx_b; }

            // get binary range
            zono_float get_bin_low() const { return bin_low; }
            zono_float get_bin_high() const { return bin_high; }

        protected:
            bool contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                                const Eigen::Vector<zono_float, -1>& b, const int iter,
                                const std::set<int>& constraints) override;

        private:
            /// binary variable indices
            std::pair<int, int> idx_b;
            zono_float bin_low = zero, bin_high = one;
        };

    } // namespace detail

    // forward declarations

    /**
     * @brief Computes the Cartesian product of two boxes Z1 and Z2.
     *
     * @param Z1 box
     * @param Z2 box
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::cartesian_product
     *
     * @note Unlike the HybZono overload, Box::operator* does NOT compute the Cartesian product.
     */
    Box cartesian_product(const Box& Z1, const Box& Z2);

    /**
     * @brief Computes the Minkowski sum of two boxes Z1 and Z2.
     *
     * @param Z1 box
     * @param Z2 box
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::minkowski_sum
     * @throws std::invalid_argument if Z1 and Z2 have different dimensions.
     */
    Box minkowski_sum(const Box& Z1, const Box& Z2);

    /**
     * @brief Computes the Pontryagin difference Z1 - Z2 of two boxes.
     *
     * @param Z1 minuend
     * @param Z2 subtrahend
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::pontry_diff
     * @throws std::invalid_argument if Z1 and Z2 have different dimensions.
     */
    Box pontry_diff(const Box& Z1, const Box& Z2);

    /**
     * @brief Projects box Z onto the dimensions specified in dims.
     *
     * @param Z box
     * @param dims vector of dimensions
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::project_onto_dims
     * @throws std::invalid_argument if any entry in dims is not a valid dimension of Z.
     */
    Box project_onto_dims(const Box& Z, const std::vector<int>& dims);

    /**
     * @brief Returns affine map R*Z + s of box Z
     *
     * @param Z box
     * @param R affine map matrix
     * @param s vector offset
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::affine_map
     * @throws std::invalid_argument if R, s, and Z have inconsistent dimensions.
     */
    Box affine_map(const Box& Z, const Eigen::SparseMatrix<zono_float>& R,
                   const Eigen::Vector<zono_float, -1>& s = Eigen::Vector<zono_float, -1>());

    /**
     * @brief Computes the interval hull of two boxes.
     *
     * @param Z1 box
     * @param Z2 box
     * @return smallest box containing both Z1 and Z2
     * @ingroup ZonoOpt_SetOperations
     * @see Box::interval_hull
     * @throws std::invalid_argument if Z1 and Z2 have different dimensions.
     */
    Box interval_hull(const Box& Z1, const Box& Z2);

    /**
     * @brief Computes the interval hull of several boxes
     *
     * @param boxes boxes for which the interval hull is to be computed
     * @return smallest box containing all boxes
     * @ingroup ZonoOpt_SetOperations
     *
     * @throws std::invalid_argument if boxes is empty or if its members have inconsistent dimensions.
     */
    Box interval_hull(const std::vector<Box>& boxes);

    /**
     * @brief Computes the generalized intersection of boxes Z1 and Z2 over the matrix R.
     *
     * @param Z1 box
     * @param Z2 box
     * @param R affine map matrix; if not provided, the identity map is used (R must then be square)
     * @param contractor_iter number of interval contractor iterations to run when R is not the identity;
     *                         ignored when R is left at its default value
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::intersect
     * @see Box::intersection_over_dims
     * @throws std::invalid_argument if R is provided and its dimensions are inconsistent with Z1 and Z2,
     * or if R is omitted and Z1, Z2 have different dimensions.
     *
     * @note An interval contractor is used to compute the over-approximation of the generalized intersection
     * when an exact intersection is not available.
     */
    Box intersection(const Box& Z1, const Box& Z2,
                     const Eigen::SparseMatrix<zono_float>& R = Eigen::SparseMatrix<zono_float>(),
                     int contractor_iter = 10);

    /**
     * @brief Computes the intersection of boxes Z1 and Z2 over the specified dimensions.
     *
     * @param Z1 box
     * @param Z2 box, must have size equal to dims.size()
     * @param dims dimensions of Z1 over which the intersection is computed
     * @return box
     * @ingroup ZonoOpt_SetOperations
     * @see Box::intersection_over_dims
     * @throws std::invalid_argument if dims.size() does not equal Z2.size(), or if any entry in dims is
     * not a valid dimension of Z1.
     */
    Box intersection_over_dims(const Box& Z1, const Box& Z2, const std::vector<int>& dims);
}