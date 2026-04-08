#ifndef ZONOOPT_BOX_HPP_
#define ZONOOPT_BOX_HPP_

/**
 * @file Box.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Box and MI_Box classes
 * @version 1.0
 * @date 2026-03-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "Interval.hpp"

namespace ZonoOpt
{
    using namespace detail;

    /**
     * @brief Box (i.e., interval vector) class
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
         */
        Box intersect(const Box& other) const;

        /**
         * @brief Box interval hull
         *
         * @param other other box
         * @return interval hull of this and other
         */
        Box interval_hull(const Box& other) const;

        /**
         * @brief Check vector continment
         * 
         * @param v vector
         * @return true if v is contained in box, false otherwise
         */
        bool contains(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Check set containment
         * 
         * @param other box to check if subset
         * @return true if other is a subset of this
         */
        bool contains_set(const Box& other) const;

        // operator overloading

        /**
         * @brief Elementwise addition
         * @param other rhs box
         * @return enclosure of this + other (elementwise)
         */
        Box operator+(const Box& other) const;

        /**
         * @brief Elementwise addition in-place
         * @param other other box
         */
        void operator+=(const Box& other);

        /**
         * @brief Elementwise addition with vector
         * @param v vector to add
         * @param box box to add
         * @return enclosure of v + box (elementwise)
         */
        Box operator+(const Eigen::Vector<zono_float, -1>& v) const;


        /**
         * @brief Elementwise addition with vector in-place
         * @param v vector
         */
        void operator+=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise addition with vector
         * @param v vector to add
         * @param box box to add
         * @return enclosure of v + box (elementwise)
         */
        friend Box operator+(const Eigen::Vector<zono_float, -1>& v, const Box& box);

        /**
         * @brief Elementwise subtraction
         * @param other rhs box
         * @return enclosure of this - other (elementwise)
         */
        Box operator-(const Box& other) const;

        /**
         * @brief Elementwise subtraction in-place
         * @param other other box
         */
        void operator-=(const Box& other);

        /**
         * @brief Elementwise subtraction with vector
         * @param v vector to subtract
         * @return enclosure of this - v (elementwise)
         */
        Box operator-(const Eigen::Vector<zono_float, -1>& v) const;

        /**
         * @brief Elementwise subtraction with vector in-place
         * @param v vector to subtract
         */
        void operator-=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise subtraction with vector
         * @param v vector
         * @param box box to subtract
         * @return enclosure of v - box (elementwise)
         */
        friend Box operator-(const Eigen::Vector<zono_float, -1>& v, const Box& box);

        /**
         * @brief Elementwise multiplication
         * @param other rhs box
         * @return enclosure of this * other (elementwise)
         */
        Box operator*(const Box& other) const;

        /**
         * @brief Elementwise multiplication in-place
         * @param other other box
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
         */
        Box operator*(const Eigen::Vector<zono_float, -1>& v) const;

        /**
         * @brief Elementwise multiplication with vector in-place
         * @param v vector to multiply
         */
        void operator*=(const Eigen::Vector<zono_float, -1>& v);

        /**
         * @brief Elementwise multiplication with vector
         * @param v vector to multiply
         * @param box box to multiply
         * @return enclosure of v * box (elementwise)
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
         */
        Box operator/(const Box& other) const;

        /**
         * @brief Elementwise division in-place
         * @param other rhs box
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
         */
        Box operator&(const Box& other) const;

        /**
         * @brief Box interval hull
         *
         * @param other box
         * @return Interval hull of this and other
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
         */
        bool contract_subset(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                             const Eigen::Vector<zono_float, -1>& b, int iter,
                             const Eigen::SparseMatrix<zono_float>& A, const std::set<int>& inds,
                             int tree_search_depth);

        /**
         * @brief Linear map of box based on interval arithmetic
         * @param A map matrix (dense)
         * @return Linear mapped box
         */
        Box linear_map(const Eigen::Matrix<zono_float, -1, -1>& A) const;

        /**
         * @brief Linear map of box based on interval arithmetic
         * @param A map matrix (sparse row major)
         * @return Linear mapped box
         */
        Box linear_map(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const;

        /**
         * @brief Linear map with vector
         * @param x vector
         * @return Interval
         */
        Interval dot(const Eigen::Vector<zono_float, -1>& x) const;

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
        friend class MI_Box;
    };

    // mixed-integer box (some bounds are fixed)
    /**
     * @brief Mixed-integer box
     *
     * Extends Box class to include variables for which may only take their upper or lower bound
     * (they may not take any value on the interior).
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
         */
        MI_Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub,
               const std::pair<int, int>& idx_b, bool zero_one_form);

        /**
         * @brief Constructor for MI_Box from vector of intervals
         * @param intervals vector intervals
         * @param idx_b indices of binary variables {start index, number of binaries}
         * @param zero_one_form flag indicating whether binary variables are in {0,1} form (true) or {-1,1} form (false)
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



}


#endif