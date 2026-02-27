#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>
#include <cmath>

using namespace ZonoOpt;

#define n_dims 3

zono_float f(const Eigen::Vector<double, n_dims>& x)
{
    return 2*std::pow(std::tan(x[0]), -2) + std::cos(x[1]/x[0])/3. + std::sin(x[0] + std::atan(x[2]))*std::sinh(x[0]) + std::exp(std::acosh(std::abs(x[1]) + 1)) - std::acos(x[0])*std::asin(x[1])/std::log(std::pow(x[2], 2));
}

Interval f_int(const Box& x)
{
    if (x.size() != n_dims)
        throw std::invalid_argument("Box must have size " + std::to_string(n_dims));

    Interval x0 = x.get_element(0);
    Interval x1 = x.get_element(1);
    Interval x2 = x.get_element(2);

    return (x0.tan().pow(-2))*2. + (x1/x0).cos()/3. + (x0 + x2.arctan()).sin()*x0.sinh() + (x1.abs() + 1).arccosh().exp() - (x0.arccos()*x1.arcsin())/(x2.pow(2)).log();
}

int main()
{
    // input intervals
    const zono_float x_min = 0.1;
    const zono_float x_max = 0.2;
    Eigen::Vector<zono_float, -1> x_lb = Eigen::Vector<zono_float, n_dims>::Constant(x_min);
    Eigen::Vector<zono_float, -1> x_ub = Eigen::Vector<zono_float, n_dims>::Constant(x_max);
    const Box x (x_lb, x_ub);

    // compute interval bounds
    const Interval f_interval = f_int(x);

    // generate random points in interval
    srand(0);
    const int n_samples = 10000;

    for (int i=0; i < n_samples; ++i)
    {
        Eigen::Vector<double, n_dims> x_sample = Eigen::Vector<double, n_dims>::Random()*(x_max - x_min)/2. + Eigen::Vector<double, n_dims>::Constant((x_max + x_min)/2.);
        const zono_float f_sample = f(x_sample);

        std::stringstream err_ss;
        err_ss << "f(" << x_sample.transpose() << ") = " << f_sample << " not in " << f_interval;
        test_assert(f_interval.contains(f_sample), err_ss.str());
    }

    return 0;
}