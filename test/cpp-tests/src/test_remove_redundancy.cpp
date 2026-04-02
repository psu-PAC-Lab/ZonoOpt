#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <iostream>

using namespace ZonoOpt;

int main() {

    // constrained zonotope that can be simplified
    Eigen::Matrix<zono_float, 1, 4> G;
    G << 3., 1., 0., 0.;
    Eigen::Vector<zono_float, 1> c;
    c << 8.;
    Eigen::Matrix<zono_float, 2, 4> A;
    A << 0.5, 0., 1., 0.,
         0., 0.5, 0., 0.5;
    Eigen::Vector<zono_float, 2> b;
    b << -0.5, -1.;

    ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // simplify
    Z.remove_redundancy();
    std::cout << Z << std::endl;

    // check that support is as expected
    std::stringstream ss;
    Eigen::Vector<zono_float, 1> d;
    d << 1.;

    zono_float sup = Z.support(d);
    ss << "Expected support = 10., got support = " << sup;
    test_assert(std::abs(sup - 10.) < zono_eps, ss.str());

    d(0) = -1.;
    sup = Z.support(d);
    ss << "Expected support = 6., got support = " << sup;
    test_assert(std::abs(sup - 6.) < zono_eps, ss.str());

    return 0;
}