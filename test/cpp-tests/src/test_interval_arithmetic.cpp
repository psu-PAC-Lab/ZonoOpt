#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>
#include <cmath>

using namespace ZonoOpt;

#define n_dims 8

zono_float f(const Eigen::Vector<double, n_dims>& x)
{
    return 2.*std::pow(std::tan(x[0]), 2) + std::cos(x[1])/3. + std::sin(x[0] + std::atan(x[2]))*std::sinh(x[7]) + std::exp(x[4]/x[3]) - std::acos(x[5])*std::asin(x[6])/std::log(x[0]);
}

Interval f_int(const Box& x)
{
    if (x.size() != n_dims)
        throw std::invalid_argument("Box must have size " + std::to_string(n_dims));

    return (x[0].tan().pow(2))*2. + x[1].cos()/3. + (x[0] + x[2].arctan()).sin()*x[7].sinh() + (x[4]/x[3]).exp() - (x[5].arccos() * x[6].arcsin())/x[0].log();
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