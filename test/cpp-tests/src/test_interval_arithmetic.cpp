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

    return 2.*(x0.tan().pow(-2)) + (x1/x0).cos()/3. + (x0 + x2.arctan()).sin()*x0.sinh() + (1. + x1.abs()).arccosh().exp() - (x0.arccos()*x1.arcsin())/(x2.pow(2)).log();
}

void test_exponent() 
{
    Interval a (0.5, 3.);
    Interval b = a.pow(456./123.);

    test_assert(!b.is_empty(), "test_exponent did not succeed");
    test_assert(std::abs(b.lower() - std::pow(0.5, 456./123.)) < 1e-6, "test_exponent lower bound incorrect");
    test_assert(std::abs(b.upper() - std::pow(3., 456./123.)) < 1e-6, "test_exponent upper bound incorrect");

    a = Interval(-3., -0.5);
    try 
    {
        b = a.pow(456./123.);
        std::cerr << "test_exponent: expected fractional power of negative interval to throw" << std::endl;
        std::exit(1);
    }
    catch (std::domain_error& e) 
    {
        // expected behavior
    }
}

void run_interval_test(zono_float x_min, zono_float x_max)
{
    // create box for input intervals
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
}

int main()
{
    // Case 1: positive range
    run_interval_test(0.1, 0.2);
    
    // Case 2: negative range, approaching 0
    run_interval_test(-0.2, -0.001);

    // Case 3: spanning 0
    run_interval_test(-1., 1.);

    // test fractional exponent
    test_exponent();

    return 0;
}