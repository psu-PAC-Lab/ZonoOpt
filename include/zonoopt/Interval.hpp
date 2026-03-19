#ifndef ZONOOPT_INTERVAL_HPP_
#define ZONOOPT_INTERVAL_HPP_

/**
 * @file Interval.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Interval class
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
        Interval* clone() const;

        // get methods

        /**
         * @brief Get lower bound
         * 
         * @return lower bound 
         */
        zono_float lower() const;

        /**
         * @brief Get upper bound
         * 
         * @return upper bound
         */
        zono_float upper() const;

        // operators

        /**
         * @brief Interval addition
         * @param other other interval
         * @return enclosure of this + other
         */
        Interval operator+(const Interval& other) const;

        /**
         * @brief Interval addition in-place
         * @param other other interval
         */
        void operator+=(const Interval& other);

        /**
         * @brief Interval addition with scalar
         *
         * @param alpha scalar
         * @return enclosure of this + alpha
         */
        Interval operator+(const zono_float alpha) const;

        /**
         * @brief Interval addition with scalar in-place
         *
         * @param alpha scalar
         */
        void operator+=(const zono_float alpha);

        /**
         * @brief Interval addition with scalar
         * 
         * @param alpha scalar
         * @param interval interval
         * @return enclosure of alpha + interval
         */
        friend Interval operator+(const zono_float alpha, const Interval& interval);

        /**
         * @brief Interval subtraction
         * @param other other interval
         * @return enclosure of this - other
         */
        Interval operator-(const Interval& other) const;

        /**
         * @brief Interval subtraction in-place
         * @param other other interval
         */
        void operator-=(const Interval& other);

        /**
         * @brief Interval subtraction with scalar
         *
         * @param alpha scalar
         * @return enclosure of this - alpha
         */
        Interval operator-(const zono_float alpha) const;

        /**
         * @brief Interval subtraction with scalar in-place
         * @param alpha scalar to subtract
         */
        void operator-=(const zono_float alpha);

        /**
         * @brief Interval subtraction with scalar
         * 
         * @param alpha scalar
         * @param interval interval
         * @return enclosure of alpha - interval
         */
        friend Interval operator-( const zono_float alpha, const Interval& interval );

        /**
         * @brief Interval multiplication
         * @param other other interval
         * @return enclosure of this * other
         */
        Interval operator*(const Interval& other) const;

        /**
         * @brief Interval multiplication in-place
         * @param other other interval
         */
        void operator*=(const Interval& other);

        /**
         * @brief Interval multiplication with scalar
         *
         * @param alpha scalar
         * @return enclosure of this * alpha
         */
        Interval operator*(const zono_float alpha) const;

        /**
         * @brief Interval multiplication with scalar in-place
         * @param alpha scalar
         */
        void operator*=(const zono_float alpha);

        /**
         * @brief Interval multiplication with scalar
         * 
         * @param alpha scalar
         * @param interval interval
         * @return enclosure of alpha * interval
         */
        friend Interval operator*(const zono_float alpha, const Interval& interval);

        /**
         * @brief Interval division
         * @param other interval to divide
         * @return enclosure of this / other
         */
        Interval operator/(const Interval& other) const;

        /**
         * @brief Interval division in-place
         * @param other interval to divide
         */
        void operator/=(const Interval& other);

        /**
         * @brief Interval division with scalar
         *
         * @param alpha scalar
         * @return enclosure of this / alpha
         */
        Interval operator/(const zono_float alpha) const;

        /**
         * @brief Interval division with scalar in-place
         * @param alpha scalar
         */
        void operator/=(const zono_float alpha);

        /**
         * @brief Interval division with scalar
         * 
         * @param alpha scalar
         * @param interval interval
         * @return enclosure of alpha / interval
         */
        friend Interval operator/(const zono_float alpha, const Interval& interval);

        /**
         * @brief Unary minus: returns -1 * this
         * 
         * @return enclosure of -this
         */
        Interval operator-() const;

        /**
         * @brief Interval intersection
         * @param other other interval
         * @return intersection of this and other
         */
        Interval operator&(const Interval& other) const;

        /**
         * @brief Interval hull
         * 
         * @param other other interval
         * @return interval hull of this and other
         */
        Interval operator|(const Interval& other) const;

        /**
         * @brief Set containment operator
         * 
         * @param other 
         * @return true if this interval is a subset of other, false otherwise
         */
        bool operator<=(const Interval& other) const;

        /**
         * @brief Set containment operator
         * 
         * @param other 
         * @return true if this interval is a superset of other, false otherwise
         */
        bool operator>=(const Interval& other) const;
        
        /**
         * @brief Set equality operator
         * 
         * @param other 
         * @return true if intervals are equal, false otherwise
         */
        bool operator==(const Interval& other) const;

        /**
         * @brief Interval inverse
         * @return enclosure of inverse
         */
        Interval inv() const;

        /**
         * @brief Interval intersection
         * @param other other interval
         * @return intersection of this and other
         */
        Interval intersect(const Interval& other) const;

        /**
         * @brief Interval hull
         * 
         * @param other other interval
         * @return interval hull of this and other
         */
        Interval interval_hull(const Interval& other) const;

        /**
         * @brief Checks whether interval contains a value
         * @param x scalar value
         * @return flag indicating if interval contains x
         */
        bool contains(const zono_float x) const;

        /**
         * @brief Set containment for intervals
         * 
         * @param other interval to check if is subset
         * @return true if this interval contains other or is equal, false otherwise
         */
        bool contains_set(const Interval& other) const;

        /**
         * @brief Checks whether interval is single-valued (i.e., width is 0 within numerical tolerance)
         * @return flag indicating if interval is single-value
         */
        bool is_single_valued() const;

        /**
         * @brief Checks whether interval is empty
         * @return flag indicating if interval is empty
         */
        bool is_empty() const;

        /**
         * @brief Get center of interval
         * @return center of interval
         */
        zono_float center() const;

        /**
         * @brief Get width of interval
         * @return width of interval
         */
        zono_float width() const;

        /**
         * @brief Get radius of interval
         * @return radius of interval
         *
         * Returns interval centered at zero with width equal to the width of the original interval
         */
        Interval radius() const;

        /**
         * @brief Get absolute value of interval
         * @return enclosure of abs(this)
         */
        Interval abs() const;

        /**
         * @brief Get square root of interval
         * @return enclosure of sqrt(this)
         */
        Interval sqrt() const;

        /**
         * @brief Interval power
         *
         * @param n power
         * @return enclosure of this^n
         */
        Interval pow(const int n) const;

        /**
         * @brief Interval power with fractional exponent
         * 
         * @param f power
         * @return enclosure of this^f 
         * 
         * Calls integer power if f is an integer within numerical tolerance.
         * Calls nth_root if f is a positive rational number within numerical tolerance.
         * Otherwise throws error.
         */
        Interval pow(const zono_float f) const;

        /**
         * @brief Interval nth root
         * @param n nth root
         * @return enclosure of root_n(this)
         */
        Interval nth_root(const int n) const;

        /**
         * @brief Compute interval containing exp(x) for all x in interval
         * @return enclosure of exp(x)
         */
        Interval exp() const;

        /**
         * @brief Compute interval containing log(x) (base e) for all x in interval
         * @return enclosure of log(x)
         */
        Interval log() const;

        /**
         * @brief Compute interval containing sin(x) for all x in interval
         * @return enclosure of sin(x)
         */
        Interval sin() const;

        /**
         * @brief Compute interval containing cos(x) for all x in interval
         * @return enclosure of cos(x)
         */
        Interval cos() const;

        /**
         * @brief Compute interval containing tan(x) for all x in interval
         * @return enclosure of tan(x)
         */
        Interval tan() const;

        /**
         * @brief Compute interval containing arcsin(x) for all x in interval
         * @return enclosure of arcsin(x)
         */
        Interval arcsin() const;

        /**
         * @brief Compute interval containing arccos(x) for all x in interval
         * @return enclosure of arccos(x)
         */
        Interval arccos() const;

        /**
         * @brief Compute interval containing arctan(x) for all x in interval
         * @return enclosure of arctan(x)
         */
        Interval arctan() const;

        /**
         * @brief Compute interval containing sinh(x) for all x in interval
         * @return enclosure of sinh(x)
         */
        Interval sinh() const;

        /**
         * @brief Compute interval containing cosh(x) for all x in interval
         * @return enclosure of cosh(x)
         */
        Interval cosh() const;

        /**
         * @brief Compute interval containing tanh(x) for all x in interval
         * @return enclosure of tanh(x)
         */
        Interval tanh() const;

        /**
         * @brief Compute interval containing arcsinh(x) for all x in interval
         * @return enclosure of arcsinh(x)
         */
        Interval arcsinh() const;

        /**
         * @brief Compute interval containing arccosh(x) for all x in interval
         * @return enclosure of arccosh(x)
         */
        Interval arccosh() const;

        /**
         * @brief Compute interval containing arctanh(x) for all x in interval
         * @return enclosure of arctanh(x)
         */
        Interval arctanh() const;

        /**
         * @brief Print method for interval
         * @return string representation of interval
         */
        std::string print() const;

        /**
         * @brief Print to ostream
         * @param os
         * @param interval
         * @return ostream reference
         */
        friend std::ostream& operator<<(std::ostream& os, const Interval& interval);


    private:

        typedef boost::numeric::interval_lib::policies<
            boost::numeric::interval_lib::save_state<rounded_transc_fixed<zono_float>>,
            boost::numeric::interval_lib::checking_base<zono_float>
        > interval_policy;

        boost::numeric::interval<zono_float, interval_policy> _val;

        explicit Interval(const boost::numeric::interval<zono_float, interval_policy>& val) : _val(val) {}
    };

} // namespace ZonoOpt

#endif
