#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <iostream>
#include <array>

using namespace ZonoOpt;

static std::string g_test_data_dir;

TEST(RemoveRedundancy, Case1_ConZonoSimplifiesToZono)
{
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
    const auto Z_rr = Z.remove_redundancy();

    Eigen::Vector<zono_float, 1> d;
    d << 1.;
    zono_float sup = Z_rr->support(d);
    EXPECT_NEAR(sup, 10., 1e-3) << "case1: expected support = 10., got " << sup;

    d(0) = -1.;
    sup = Z_rr->support(d);
    EXPECT_NEAR(sup, -4., 1e-3) << "case1: expected support = -4., got " << sup;

    EXPECT_TRUE(Z_rr->is_zono()) << "case1: expected set to be a zonotope after removing redundancy.";
}

TEST(RemoveRedundancy, Case2_HybZonoSupportPreserved)
{
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

    Eigen::Vector<zono_float, 1> d;
    std::array<zono_float, 2> sup_before;
    d << 1.;
    sup_before[0] = Z.support(d);
    d(0) = -1.;
    sup_before[1] = Z.support(d);

    const auto Z_rr = Z.remove_redundancy();

    std::array<zono_float, 2> sup_after;
    d(0) = 1.;
    sup_after[0] = Z_rr->support(d);
    d(0) = -1.;
    sup_after[1] = Z_rr->support(d);

    for (int i=0; i<2; i++)
    {
        EXPECT_NEAR(sup_after[i], sup_before[i], 1e-3)
            << "case2: expected support = " << sup_before[i] << ", got " << sup_after[i];
    }
}

TEST(RemoveRedundancy, Case3_ConZonoSupportPreserved)
{
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

    Eigen::Vector<zono_float, 1> d;
    std::array<zono_float, 2> sup_before;
    d << 1.;
    sup_before[0] = Z.support(d);
    d(0) = -1.;
    sup_before[1] = Z.support(d);

    const auto Z_rr = Z.remove_redundancy();

    std::array<zono_float, 2> sup_after;
    d(0) = 1.;
    sup_after[0] = Z_rr->support(d);
    d(0) = -1.;
    sup_after[1] = Z_rr->support(d);

    for (int i=0; i<2; i++)
    {
        EXPECT_NEAR(sup_after[i], sup_before[i], 1e-3)
            << "case3: expected support = " << sup_before[i] << ", got " << sup_after[i];
    }
}

TEST(RemoveRedundancy, Case4_InfeasibleBecomesEmptySet)
{
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
    const auto Z_rr = Z.remove_redundancy();

    std::stringstream ss;
    ss << "Expected EmptySet, got " << *Z_rr;
    EXPECT_TRUE(Z_rr->is_empty_set()) << ss.str();
}

TEST(RemoveRedundancy, NonEmptyFromFile_test5)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/remove_redundancy/Z_test5.json";

    const ZonoPtr Z = from_json(filename);
    OptSettings settings;
    settings.verbose = true;

    const auto Z_rr = Z->remove_redundancy();
    std::cout << *Z_rr << std::endl;

    std::stringstream ss;
    ss << "Expected non-empty set, got " << *Z_rr;
    EXPECT_FALSE(Z_rr->is_empty(settings)) << ss.str();
}

TEST(RemoveRedundancy, NonEmptyFromFile_test6)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/remove_redundancy/Z_test6.json";

    const ZonoPtr Z = from_json(filename);
    OptSettings settings;
    settings.verbose = true;

    const auto Z_rr = Z->remove_redundancy();
    std::cout << *Z_rr << std::endl;

    std::stringstream ss;
    ss << "Expected non-empty set, got " << *Z_rr;
    EXPECT_FALSE(Z_rr->is_empty(settings)) << ss.str();
}

TEST(RemoveRedundancy, NonEmptyFromFile_test7)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/remove_redundancy/Z_test7.json";

    const ZonoPtr Z = from_json(filename);
    OptSettings settings;
    settings.verbose = true;

    const auto Z_rr = Z->remove_redundancy();
    std::cout << *Z_rr << std::endl;

    std::stringstream ss;
    ss << "Expected non-empty set, got " << *Z_rr;
    EXPECT_FALSE(Z_rr->is_empty(settings)) << ss.str();
}

TEST(RemoveRedundancy, RandomConZono100)
{
    std::mt19937 rand_gen(0);

    OptSettings settings;
    settings.eps_prim = 1e-3;
    settings.eps_dual = 1e-3;
    settings.rho = 1.;

    for (int i=0; i<100; ++i)
    {
        ConZono Z = random_conzono(2, 30, 10, 0.1, 0., 1., rand_gen);

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
        catch (std::runtime_error&) { continue; }
        catch (std::invalid_argument&) { continue; }
        catch (std::exception& err)
        {
            ADD_FAILURE() << err.what();
            continue;
        }

        std::uniform_real_distribution<double> form_dist(0., 1.);
        if (form_dist(rand_gen) < 0.5)
            Z.convert_form();

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

        for (int j=0; j<4; ++j)
        {
            const bool close = std::abs(sup_before[j] - sup_after[j])/std::abs(sup_before[j]) < 1e-2
                            || std::abs(sup_before[j] - sup_after[j]) < 1e-1;
            EXPECT_TRUE(close)
                << "Random ConZono iter " << i << ": expected support = " << sup_before[j]
                << ", got support = " << sup_after[j]
                << "\n  Z before: " << Z << "\n  Z after: " << *Z_rr;
        }
    }
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc >= 2) g_test_data_dir = argv[1];
    return RUN_ALL_TESTS();
}
