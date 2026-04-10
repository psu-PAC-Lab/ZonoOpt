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

    const ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // simplify
    const auto Z_rr = Z.remove_redundancy();

    // check that support is as expected
    std::stringstream ss;
    Eigen::Vector<zono_float, 1> d;
    d << 1.;

    zono_float sup = Z_rr->support(d);
    ss << "case1: expected support = 10., got support = " << sup << std::endl;
    test_assert(std::abs(sup - 10.) < 1e-3, ss.str());
    ss.str("");

    d(0) = -1.;
    sup = Z_rr->support(d);
    ss << "case1: expected support = -4., got support = " << sup << std::endl;
    test_assert(std::abs(sup - -4.) < 1e-3, ss.str());

    // set should be a zonotope
    test_assert(Z_rr->is_zono(), "case1: expected set to be a zonotope after removing redundancy.");
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

    HybZono Z (Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);

    // get support
    std::array<zono_float, 2> sup_before;

    Eigen::Vector<zono_float, 1> d;
    d << 1.;
    sup_before[0] = Z.support(d);

    d(0) = -1.;
    sup_before[1] = Z.support(d);

    // simplify
    const auto Z_rr = Z.remove_redundancy();

    // check that support is as expected
    std::stringstream ss;
    std::array<zono_float, 2> sup_after;
    d(0) = 1.;
    sup_after[0] = Z_rr->support(d);
    d(0) = -1.;
    sup_after[1] = Z_rr->support(d);

    for (int i=0; i<2; i++)
    {
        ss << "case2: expected support = " << sup_before[i] << ", got support = " << sup_after[i] << std::endl;
        test_assert(std::abs(sup_before[i] - sup_after[i]) < 1e-3, ss.str());
        ss.str("");
    }
}

void test3()
{
    // constrained zonotope
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
    std::array<zono_float, 2> sup_before;

    Eigen::Vector<zono_float, 1> d;
    d << 1.;
    sup_before[0] = Z.support(d);

    d(0) = -1.;
    sup_before[1] = Z.support(d);

    // simplify
    const auto Z_rr = Z.remove_redundancy();

    // check that support is as expected
    std::stringstream ss;
    std::array<zono_float, 2> sup_after;
    d(0) = 1.;
    sup_after[0] = Z_rr->support(d);
    d(0) = -1.;
    sup_after[1] = Z_rr->support(d);

    for (int i=0; i<2; i++)
    {
        ss << "case3: expected support = " << sup_before[i] << ", got support = " << sup_after[i] << std::endl;
        test_assert(std::abs(sup_before[i] - sup_after[i]) < 1e-3, ss.str());
        ss.str("");
    }
}

void test4()
{
    // infeasible constrained zonotope
    Eigen::Matrix<zono_float, 2, 2> G;
    G << 1., 1.,
         1., 2.;
    Eigen::Vector<zono_float, 2> c;
    c << 1., 2.;
    Eigen::Matrix<zono_float, 1, 2> A;
    A << 1., 1.;
    Eigen::Vector<zono_float, 1> b;
    b << 3.;

    const ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // remove redundancy
    const auto Z_rr = Z.remove_redundancy();

    // check that result is EmptySet object
    std::stringstream ss;
    ss << "Expected EmptySet, got " << *Z_rr << std::endl;
    test_assert(Z_rr->is_empty_set(), ss.str());
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
    catch (std::runtime_error& err)
    {
        return;
    }
    catch (std::invalid_argument& err)
    {
        return;
    }
    catch (std::exception& err)
    {
        std::cerr << err.what() << std::endl;
        std::exit(1);
    }

    // randomly convert form
    std::uniform_real_distribution<double> form_dist(0., 1.);
    if (form_dist(rand_gen) < 0.5)
        Z.convert_form();

    // get support after simplifying
    const auto Z_rr = Z.remove_redundancy();
    std::array<zono_float, 4> sup_after;

    d << 1., 0.;
    sup_after[0] = Z_rr->support(d, settings);
    d << -1., 0.;
    sup_after[1] = Z_rr->support(d, settings);
    d << 0., 1.;
    sup_after[2] = Z_rr->support(d, settings);
    d << 0., -1.;
    sup_after[3] = Z_rr->support(d, settings);

    // make sure all close
    std::stringstream ss;
    for (int i=0; i<4; ++i)
    {
        ss << "Random ConZono: expected support = " << sup_before[i] << ", got support = " << sup_after[i] << std::endl;
        ss << "  Z before simplifying: " << Z << std::endl;
        ss << "  Z after simplifying: " << *Z_rr << std::endl;

        test_assert(std::abs(sup_before[i] - sup_after[i])/std::abs(sup_before[i]) < 1e-2 || std::abs(sup_before[i] - sup_after[i]) < 1e-1, ss.str());
        ss.str("");
    }
}

int main()
{
    // manually specified test
    test1();
    test2();
    test3();
    test4();

    // random constrained and hybrid zonotopes
    std::mt19937 rand_gen(0);

    for (int i=0; i<100; ++i)
    {
        test_random_conzono(rand_gen);
    }

    return 0;
}