#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <iostream>
#include <array>

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
    ss.str("");

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
    ss.str("");

    // check that support is as expected
    d(0) = 1.;
    const zono_float sup1_2 = Zh.support(d);
    ss << "case3: expected support = " << sup1_1 << ", got support = " << sup1_2 << std::endl;
    test_assert(std::abs(sup1_1 - sup1_2) < 1e-3, ss.str());
    ss.str("");

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
    ss.str("");

    d(0) = -1.;
    const zono_float sup2_2 = Z.support(d);
    ss << "case3: expected support = " << sup2_1 << ", got support = " << sup2_2 << std::endl;
    test_assert(std::abs(sup2_1 - sup2_2) < 1e-3, ss.str());
}

void test_random_conzono(std::mt19937& rand_gen)
{
    ConZono Z = random_conzono(2, 30, 10, 0.1, 0., 1., rand_gen);

    // get support before simplifying
    OptSettings settings;
    settings.eps_prim = 1e-3;
    settings.eps_dual = 1e-3;
    settings.rho = 1.;
    std::array<zono_float, 4> sup_before;

    Eigen::Vector<zono_float, 2> d;
    try
    {
        d << 1., 0.;
        sup_before[0] = Z.support(d, settings);
        d << -1., 0.;
        sup_before[1] = Z.support(d, settings);
        d << 0., 1.;
        sup_before[2] = Z.support(d, settings);
        d << 0., -1.;
        sup_before[3] = Z.support(d, settings);
    }
    catch (std::exception& err)
    {
        return;
    }

    // randomly convert form
    // std::uniform_real_distribution<double> form_dist(0., 1.);
    // if (form_dist(rand_gen) < 0.5)
    //     Z.convert_form();

    std::string Z_before_str = Z.print();

    // get support after simplifying
    Z.remove_redundancy();
    std::array<zono_float, 4> sup_after;

    d << 1., 0.;
    sup_after[0] = Z.support(d, settings);
    d << -1., 0.;
    sup_after[1] = Z.support(d, settings);
    d << 0., 1.;
    sup_after[2] = Z.support(d, settings);
    d << 0., -1.;
    sup_after[3] = Z.support(d, settings);

    // make sure all close
    std::stringstream ss;
    for (int i=0; i<4; ++i)
    {
        ss << "Random ConZono: expected support = " << sup_before[i] << ", got support = " << sup_after[i] << std::endl;
        ss << "  Z before simplifying: " << Z_before_str << std::endl;
        ss << "  Z after simplifying: " << Z << std::endl;

        test_assert(std::abs(sup_before[i] - sup_after[i])/std::abs(sup_before[i]) < 1e-2 || std::abs(sup_before[i] - sup_after[i]) < 1e-3, ss.str());
        ss.str("");
    }
}

int main()
{
    // manually specified test
    test1();
    test2();
    test3();

    // random constrained zonotopes
    std::mt19937 rand_gen(0);

    for (int i=0; i<500; ++i)
    {
        if (i == 301)
            std::cout << "DEBUG" << std::endl;

        test_random_conzono(rand_gen);
    }


    return 0;
}