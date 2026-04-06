#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <iostream>

using namespace ZonoOpt;

void test1()
{
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

    // check that support is as expected
    std::stringstream ss;
    Eigen::Vector<zono_float, 1> d;
    d << 1.;

    zono_float sup = Z.support(d);
    ss << "case1: expected support = 10., got support = " << sup << std::endl;
    test_assert(std::abs(sup - 10.) < 1e-3, ss.str());
    ss.clear();

    d(0) = -1.;
    sup = Z.support(d);
    ss << "case1: expected support = -4., got support = " << sup << std::endl;
    test_assert(std::abs(sup - -4.) < 1e-3, ss.str());
}

void test2()
{
    // hybrid zonotope
    Eigen::Matrix<zono_float, 1, 2> Gc;
    Gc << 3., 1.;
    Eigen::Matrix<zono_float, 1, 2> Gb;
    Gb << 0., 0.;
    Eigen::Vector<zono_float, 1> c;
    c << 8.;
    Eigen::Matrix<zono_float, 2, 2> Ac;
    Ac << 0.5, 0.,
         0., 0.5;
    Eigen::Matrix<zono_float, 2, 2> Ab;
    Ab << 1., 0.,
          0., 0.5;
    Eigen::Vector<zono_float, 2> b;
    b << -0.5, -1.;

    HybZono Zh (Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);

    // get leaves
    auto leaves_before_simplifying = Zh.get_leaves(false);

    // get support
    Eigen::Vector<zono_float, 1> d;
    d << 1.;
    const zono_float sup1_1 = Zh.support(d);

    d(0) = -1.;
    const zono_float sup2_1 = Zh.support(d);

    // simplify
    Zh.remove_redundancy();

    // get leaves
    auto leaves_after_simplifying = Zh.get_leaves(false);

    std::stringstream ss;
    ss << "case2: expected " << leaves_before_simplifying.size() << " leaves, got " << leaves_after_simplifying.size() << std::endl;
    test_assert(leaves_before_simplifying.size() == leaves_after_simplifying.size(), ss.str());
    ss.clear();

    // check that support is as expected
    d(0) = 1.;
    const zono_float sup1_2 = Zh.support(d);
    ss << "case3: expected support = " << sup1_1 << ", got support = " << sup1_2 << std::endl;
    test_assert(std::abs(sup1_1 - sup1_2) < 1e-3, ss.str());
    ss.clear();

    d(0) = -1.;
    const zono_float sup2_2 = Zh.support(d);
    ss << "case3: expected support = " << sup2_1 << ", got support = " << sup2_2 << std::endl;
    test_assert(std::abs(sup2_1 - sup2_2) < 1e-3, ss.str());
}

void test3()
{
    // constrained zonotope that can be simplified
    Eigen::Matrix<zono_float, 1, 4> G;
    G << 3., 1., 0., 0.;
    Eigen::Vector<zono_float, 1> c;
    c << 8.;
    Eigen::Matrix<zono_float, 2, 4> A;
    A << 0.5, 0.1, 1., 0.,
         0., 0.5, 0., 0.5;
    Eigen::Vector<zono_float, 2> b;
    b << -0.5, -1.;

    ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // get support
    Eigen::Vector<zono_float, 1> d;
    d << 1.;
    const zono_float sup1_1 = Z.support(d);

    d(0) = -1.;
    const zono_float sup2_1 = Z.support(d);

    // simplify
    Z.remove_redundancy();

    // check that support is as expected
    std::stringstream ss;

    d(0) = 1.;
    const zono_float sup1_2 = Z.support(d);
    ss << "case3: expected support = " << sup1_1 << ", got support = " << sup1_2 << std::endl;
    test_assert(std::abs(sup1_1 - sup1_2) < 1e-3, ss.str());
    ss.clear();

    d(0) = -1.;
    const zono_float sup2_2 = Z.support(d);
    ss << "case3: expected support = " << sup2_1 << ", got support = " << sup2_2 << std::endl;
    test_assert(std::abs(sup2_1 - sup2_2) < 1e-3, ss.str());
}

int main()
{
    test1();
    test2();
    test3();


    return 0;
}