#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

        Interval* Interval::clone() const
        {
            return new Interval(*this);
        }

        zono_float Interval::lower() const { return _val.lower(); }

        zono_float Interval::upper() const { return _val.upper(); }

        // operators

        Interval Interval::operator+(const Interval& other) const
        {
            return Interval(this->_val + other._val);
        }

        void Interval::operator+=(const Interval& other)
        {
            this->_val += other._val;
        }

        Interval Interval::operator+(const zono_float alpha) const
        {
            return Interval(this->_val + alpha);
        }

        void Interval::operator+=(const zono_float alpha)
        {
            this->_val += alpha;
        }

        Interval operator+(const zono_float alpha, const Interval& interval)
        {
            return interval + alpha;
        }

        Interval Interval::operator-(const Interval& other) const
        {
            return Interval(this->_val - other._val);
        }

        void Interval::operator-=(const Interval& other)
        {
            this->_val -= other._val;
        }

        Interval Interval::operator-(const zono_float alpha) const
        {
            return Interval(this->_val - alpha);
        }

        void Interval::operator-=(const zono_float alpha)
        {
            this->_val -= alpha;
        }

        Interval operator-( const zono_float alpha, const Interval& interval )
        {
            return Interval(alpha - interval._val);
        }

        Interval Interval::operator*(const Interval& other) const
        {
            return Interval(this->_val * other._val);
        }

        void Interval::operator*=(const Interval& other)
        {
            this->_val *= other._val;
        }

        Interval Interval::operator*(const zono_float alpha) const
        {
            return Interval(this->_val * alpha);
        }

        void Interval::operator*=(const zono_float alpha)
        {
            this->_val *= alpha;
        }

        Interval operator*(const zono_float alpha, const Interval& interval)
        {
            return interval*alpha;
        }

        Interval Interval::operator/(const Interval& other) const
        {
            return Interval(this->_val / other._val);
        }

        void Interval::operator/=(const Interval& other)
        {
            this->_val /= other._val;
        }

        Interval Interval::operator/(const zono_float alpha) const
        {
            return Interval(this->_val / alpha);
        }

        void Interval::operator/=(const zono_float alpha)
        {
            this->_val /= alpha;
        }

        Interval operator/(const zono_float alpha, const Interval& interval)
        {
            return Interval(alpha / interval._val);
        }

        Interval Interval::operator-() const
        {
            return Interval(-this->_val);
        }

        Interval Interval::operator&(const Interval& other) const
        {
            return this->intersect(other);
        }

        Interval Interval::operator|(const Interval& other) const
        {
            return this->interval_hull(other);
        }

        bool Interval::operator<=(const Interval& other) const
        {
            return other.contains_set(*this);
        }

        bool Interval::operator>=(const Interval& other) const
        {
            return this->contains_set(other);
        }

        bool Interval::operator==(const Interval& other) const
        {
            return this->contains_set(other) && other.contains_set(*this);
        }

        Interval Interval::inv() const
        {
            return Interval(boost::numeric::interval_lib::multiplicative_inverse(_val));
        }

        Interval Interval::intersect(const Interval& other) const
        {
            return Interval(boost::numeric::intersect(_val, other._val));
        }

        Interval Interval::interval_hull(const Interval& other) const
        {
            return Interval(boost::numeric::hull(_val, other._val));
        }

        bool Interval::contains(const zono_float x) const
        {
            return x >= this->lower() - zono_eps && x <= this->upper() + zono_eps;
        }

        bool Interval::contains_set(const Interval& other) const
        {
            return this->contains(other.lower()) && this->contains(other.upper());
        }

        bool Interval::is_single_valued() const
        {
            return this->width() < zono_eps;
        }

        bool Interval::is_empty() const
        {
            return boost::numeric::empty(_val);
        }

        zono_float Interval::center() const
        {
            return boost::numeric::median(_val);
        }

        zono_float Interval::width() const
        {
            return _val.upper() - _val.lower();
        }

        Interval Interval::radius() const
        {
            const zono_float r = this->width() / two;
            return Interval(-r, r);
        }

        Interval Interval::abs() const
        {
            return Interval(boost::numeric::abs(_val));
        }

        Interval Interval::sqrt() const
        {
            return Interval(boost::numeric::sqrt(_val));
        }

        Interval Interval::pow(const int n) const
        {
            return Interval(boost::numeric::pow(_val, n));
        }

        Interval Interval::pow(const zono_float f) const
        {
            // integer-valued case
            if (std::abs(f - std::round(f)) < zono_eps)
                return this->pow(static_cast<int>(std::round(f)));

            // negative intervals not supported
            if (this->lower() < zero)
                throw std::domain_error("Negative intervals not supported for non-integer powers.");

            // get rational exponent
            const auto [numerator, denominator] = get_rational(static_cast<double>(f));
            return this->pow(numerator).nth_root(denominator);
        }

        Interval Interval::nth_root(const int n) const
        {
            return Interval(boost::numeric::nth_root(_val, n));
        }

        Interval Interval::exp() const
        {
            return Interval(boost::numeric::exp(_val));
        }

        Interval Interval::log() const
        {
            return Interval(boost::numeric::log(_val));
        }

        Interval Interval::sin() const
        {
            return Interval(boost::numeric::sin(_val));
        }

        Interval Interval::cos() const
        {
            return Interval(boost::numeric::cos(_val));
        }

        Interval Interval::tan() const
        {
            return Interval(boost::numeric::tan(_val));
        }

        Interval Interval::arcsin() const
        {
            return Interval(boost::numeric::asin(_val));
        }

        Interval Interval::arccos() const
        {
            return Interval(boost::numeric::acos(_val));
        }

        Interval Interval::arctan() const
        {
            return Interval(boost::numeric::atan(_val));
        }

        Interval Interval::sinh() const
        {
            return Interval(boost::numeric::sinh(_val));
        }

        Interval Interval::cosh() const
        {
            return Interval(boost::numeric::cosh(_val));
        }

        Interval Interval::tanh() const
        {
            return Interval(boost::numeric::tanh(_val));
        }

        Interval Interval::arcsinh() const
        {
            return Interval(boost::numeric::asinh(_val));
        }

        Interval Interval::arccosh() const
        {
            return Interval(boost::numeric::acosh(_val));
        }

        Interval Interval::arctanh() const
        {
            return Interval(boost::numeric::atanh(_val));
        }

        std::string Interval::print() const
        {
            return "Interval: [" + std::to_string(lower()) + ", " + std::to_string(upper()) + "]";
        }

        std::ostream& operator<<(std::ostream& os, const Interval& interval)
        {
            os << interval.print();
            return os;
        }

        std::pair<int, int> Interval::get_rational(zono_float x)
        {
            // get upper and lower bound
            int n_l = static_cast<int>(std::floor(x)); // numerator, lower
            int n_u = static_cast<int>(std::ceil(x)); // numerator, upper
            int d_l = 1; // denominator, lower
            int d_u = 1; // denominator, upper

            // get point in between
            int n = n_l + n_u;
            int d = d_l + d_u;
            zono_float n_d = static_cast<zono_float>(n) / d;

            // iterate until convergence
            while (std::abs(n_d - x) > zono_eps) 
            {
                if (x > n_d)
                {
                    n_l = n;
                    d_l = d;
                }
                else
                {
                    n_u = n;
                    d_u = d;
                }
                n = n_l + n_u;
                d = d_l + d_u;
                n_d = static_cast<zono_float>(n) / d;
            }

            return {n, d};
        }
}
