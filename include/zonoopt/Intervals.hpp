#ifndef ZONOOPT_INTERVAL_UTILITIES_HPP_
#define ZONOOPT_INTERVAL_UTILITIES_HPP_

/**
 * @file Intervals.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Intervals and associated classes.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <string>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <boost/numeric/interval.hpp>

namespace ZonoOpt
{
    namespace detail
    {
        // need this to get inverse hyperbolics to compile with MSVC
        template<class T>
        struct rounded_transc_fixed : boost::numeric::interval_lib::rounded_transc_std<T> {
            static T asinh_down(const T& x) { using std::asinh; return asinh(x); }
            static T asinh_up  (const T& x) { using std::asinh; return asinh(x); }
            static T acosh_down(const T& x) { using std::acosh; return acosh(x); }
            static T acosh_up  (const T& x) { using std::acosh; return acosh(x); }
            static T atanh_down(const T& x) { using std::atanh; return atanh(x); }
            static T atanh_up  (const T& x) { using std::atanh; return atanh(x); }
        };
    }

    using namespace detail;


    /**
     * @brief Interval class
     * 
     * Wraps boost::numeric::interval
     */
    class Interval
    {
    public:
        // constructor

        /**
         * @brief Default constructor
         */
        Interval() : _val(zero, zero) {}

        /**
         * @brief Interval constructor
         * @param y_min lower bound
         * @param y_max upper bound
         */
        Interval(const zono_float y_min, const zono_float y_max) : _val(y_min, y_max) {}

        /**
         * @brief Clone Interval object
         * @return clone of object
         */
        Interval* clone() const
        {
            return new Interval(*this);
        }

        // get methods

        /**
         * @brief Get lower bound
         * 
         * @return lower bound 
         */
        zono_float lb() const { return _val.lower(); }

        /**
         * @brief Get upper bound
         * 
         * @return upper bound
         */
        zono_float ub() const { return _val.upper(); }

        // operators

        /**
         * @brief Interval addition
         * @param other other interval
         * @return enclosure of this + other
         */
        Interval operator+(const Interval& other) const
        {
            return Interval(this->_val + other._val);
        }

        /**
         * @brief Interval addition with scalar
         *
         * @param alpha scalar
         * @return enclosure of this + alpha
         */
        Interval operator+(const zono_float alpha) const
        {
            return Interval(this->_val + alpha);
        }

        /**
         * @brief Interval subtraction
         * @param other other interval
         * @return enclosure of this - other
         */
        Interval operator-(const Interval& other) const
        {
            return Interval(this->_val - other._val);
        }

        /**
         * @brief Interval subtraction with scalar
         *
         * @param alpha scalar
         * @return enclosure of this - alpha
         */
        Interval operator-(const zono_float alpha) const
        {
            return Interval(this->_val - alpha);
        }

        /**
         * @brief Interval multiplication
         * @param other other interval
         * @return enclosure of this * other
         */
        Interval operator*(const Interval& other) const
        {
            return Interval(this->_val * other._val);
        }

        /**
         * @brief Interval multiplication with scalar
         *
         * @param alpha scalar
         * @return enclosure of this * alpha
         */
        Interval operator*(const zono_float alpha) const
        {
            return Interval(this->_val * alpha);
        }

        /**
         * @brief Interval division
         * @param other other interval
         * @return enclosure of this / other
         */
        Interval operator/(const Interval& other) const
        {
            return Interval(this->_val / other._val);
        }

        /**
         * @brief Interval division with scalar
         *
         * @param alpha scalar
         * @return enclosure of this / alpha
         */
        Interval operator/(const zono_float alpha) const
        {
            return Interval(this->_val / alpha);
        }

        /**
         * @brief Interval inverse
         * @return enclosure of inverse
         */
        Interval inv() const
        {
            return Interval(boost::numeric::interval_lib::multiplicative_inverse(_val));
        }

        /**
         * @brief Interval intersection
         * @param other other interval
         * @return intersection of this and other
         */
        Interval intersect(const Interval& other) const
        {
            return Interval(boost::numeric::intersect(_val, other._val));
        }

        /**
         * @brief Checks whether interval contains a value
         * @param x scalar value
         * @return flag indicating if interval contains x
         */
        bool contains(const zono_float x) const
        {
            return x >= this->lb() - zono_eps && x <= this->ub() + zono_eps;
        }

        /**
         * @brief Checks whether interval is single-valued (i.e., width is 0 within numerical tolerance)
         * @return flag indicating if interval is single-value
         */
        bool is_single_valued() const
        {
            return this->width() < zono_eps;
        }

        /**
         * @brief Checks whether interval is empty
         * @return flag indicating if interval is empty
         */
        bool is_empty() const
        {
            return boost::numeric::empty(_val);
        }

        /**
         * @brief Get center of interval
         * @return center of interval
         */
        zono_float center() const
        {
            return boost::numeric::median(_val);
        }

        /**
         * @brief Get width of interval
         * @return width of interval
         */
        zono_float width() const
        {
            return _val.upper() - _val.lower();
        }

        /**
         * @brief Get radius of interval
         * @return radius of interval
         *
         * Returns interval centered at zero with width equal to the width of the original interval
         */
        Interval radius() const
        {
            const zono_float r = this->width() / two;
            return Interval(-r, r);
        }


        /**
         * @brief Get absolute value of interval
         * @return enclosure of abs(this)
         */
        Interval abs() const
        {
            return Interval(boost::numeric::abs(_val));
        }

        /**
         * @brief Get square root of interval
         * @return enclosure of sqrt(this)
         */
        Interval sqrt() const
        {
            return Interval(boost::numeric::sqrt(_val));
        }

        /**
         * @brief Interval power
         *
         * @param n power
         * @return enclosure of this^n
         */
        Interval pow(const int n) const
        {
            return Interval(boost::numeric::pow(_val, n));
        }

        /**
         * @brief Interval nth root
         * @param n nth root
         * @return enclosure of root_n(this)
         */
        Interval nth_root(const int n) const
        {
            return Interval(boost::numeric::nth_root(_val, n));
        }

        /**
         * @brief Compute interval containing exp(x) for all x in interval
         * @return enclosure of exp(x)
         */
        Interval exp() const
        {
            return Interval(boost::numeric::exp(_val));
        }

        /**
         * @brief Compute interval containing log(x) (base e) for all x in interval
         * @return enclosure of log(x)
         */
        Interval log() const
        {
            return Interval(boost::numeric::log(_val));
        }

        /**
         * @brief Compute interval containing sin(x) for all x in interval
         * @return enclosure of sin(x)
         */
        Interval sin() const
        {
            return Interval(boost::numeric::sin(_val));
        }

        /**
         * @brief Compute interval containing cos(x) for all x in interval
         * @return enclosure of cos(x)
         */
        Interval cos() const
        {
            return Interval(boost::numeric::cos(_val));
        }

        /**
         * @brief Compute interval containing tan(x) for all x in interval
         * @return enclosure of tan(x)
         */
        Interval tan() const
        {
            return Interval(boost::numeric::tan(_val));
        }

        /**
         * @brief Compute interval containing arcsin(x) for all x in interval
         * @return enclosure of arcsin(x)
         */
        Interval arcsin() const
        {
            return Interval(boost::numeric::asin(_val));
        }

        /**
         * @brief Compute interval containing arccos(x) for all x in interval
         * @return enclosure of arccos(x)
         */
        Interval arccos() const
        {
            return Interval(boost::numeric::acos(_val));
        }

        /**
         * @brief Compute interval containing arctan(x) for all x in interval
         * @return enclosure of arctan(x)
         */
        Interval arctan() const
        {
            return Interval(boost::numeric::atan(_val));
        }

        /**
         * @brief Compute interval containing sinh(x) for all x in interval
         * @return enclosure of sinh(x)
         */
        Interval sinh() const
        {
            return Interval(boost::numeric::sinh(_val));
        }

        /**
         * @brief Compute interval containing cosh(x) for all x in interval
         * @return enclosure of cosh(x)
         */
        Interval cosh() const
        {
            return Interval(boost::numeric::cosh(_val));
        }

        /**
         * @brief Compute interval containing tanh(x) for all x in interval
         * @return enclosure of tanh(x)
         */
        Interval tanh() const
        {
            return Interval(boost::numeric::tanh(_val));
        }

        /**
         * @brief Compute interval containing arcsinh(x) for all x in interval
         * @return enclosure of arcsinh(x)
         */
        Interval arcsinh() const
        {
            return Interval(boost::numeric::asinh(_val));
        }

        /**
         * @brief Compute interval containing arccosh(x) for all x in interval
         * @return enclosure of arccosh(x)
         */
        Interval arccosh() const
        {
            return Interval(boost::numeric::acosh(_val));
        }

        /**
         * @brief Compute interval containing arctanh(x) for all x in interval
         * @return enclosure of arctanh(x)
         */
        Interval arctanh() const
        {
            return Interval(boost::numeric::atanh(_val));
        }

        /**
         * @brief Print method for interval
         * @return string representation of interval
         */
        std::string print() const
        {
            return "Interval: [" + std::to_string(lb()) + ", " + std::to_string(ub()) + "]";
        }

        /**
         * @brief Print to ostream
         * @param os
         * @param interval
         * @return ostream reference
         */
        friend std::ostream& operator<<(std::ostream& os, const Interval& interval)
        {
            os << interval.print();
            return os;
        }


    private:

        typedef boost::numeric::interval_lib::policies<
            boost::numeric::interval_lib::save_state<rounded_transc_fixed<zono_float>>,
            boost::numeric::interval_lib::checking_base<zono_float>
        > interval_policy;

        boost::numeric::interval<zono_float, interval_policy> _val;

        explicit Interval(const boost::numeric::interval<zono_float, interval_policy>& val) : _val(val) {}
    };


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

        // operator overloading

        /**
         * @brief Elementwise addition
         * @param other rhs box
         * @return enclosure of this + other (elementwise)
         */
        Box operator+(const Box& other) const;

        /**
         * @brief Elementwise subtraction
         * @param other rhs box
         * @return enclosure of this - other (elementwise)
         */
        Box operator-(const Box& other) const;

        /**
         * @brief Elementwise multiplication
         * @param other rhs box
         * @return enclosure of this * other (elementwise)
         */
        Box operator*(const Box& other) const;

        /**
         * @brief Elementwise multiplication with scalar
         * @param alpha scalar multiplier
         * @return enclosure of alpha * this (elementwise)
         */
        Box operator*(zono_float alpha) const;

        /**
         * @brief Elementwise division
         * @param other rhs box
         * @return enclosure of this / other (elementwise)
         */
        Box operator/(const Box& other) const;

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

    protected:
        bool contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                             const Eigen::Vector<zono_float, -1>& b, const int iter,
                             const std::set<int>& constraints) override;

    private:
        /// binary variable indices
        std::pair<int, int> idx_b;
        zono_float bin_low = zero, bin_high = one;
    };


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
         * @brief IntervalMatrix constructor using sparse lower and upper bound matrices
         * @param mat_lb lower bound matrix
         * @param mat_ub upper bound matrix
         */
        IntervalMatrix(const Eigen::SparseMatrix<zono_float>& mat_lb,
                       const Eigen::SparseMatrix<zono_float>& mat_ub);

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
        IntervalMatrix(const Eigen::Matrix<Interval, -1, -1>& mat);

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
         * @brief IntervalMatrix multiplication with matrix
         * @param A rhs matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const;

        /**
         * @brief IntervalMatrix multiplication with another IntervalMatrix
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator*(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix addition
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator+(const IntervalMatrix& other) const;

        /**
         * @brief IntervalMatrix subtraction
         * @param other rhs interval matrix
         * @return resulting interval matrix
         */
        IntervalMatrix operator-(const IntervalMatrix& other) const;

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
} // namespace ZonoOpt

#endif
