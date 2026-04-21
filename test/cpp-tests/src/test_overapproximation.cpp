#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

static std::string g_test_data_dir;

TEST(Overapproximation, ReduceOrder)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/overapproximation/rand_zono.json";

    ZonoPtr Z = from_json(filename);

    Zono* Z_zono_ptr = dynamic_cast<Zono*>(Z.get());
    ASSERT_NE(Z_zono_ptr, nullptr) << "Expected set to be a zonotope for reduce order test.";
    Z.release();
    std::unique_ptr<Zono> Z_zono (Z_zono_ptr);

    Eigen::Vector<zono_float, -1> p (Z_zono->get_n());
    p.setZero();
    const auto p_proj = Z_zono->project_point(p);

    const auto Zr = Z_zono->reduce_order(10);

    EXPECT_TRUE(Zr->contains_point(p_proj)) << "Reduced zonotope does not contain projected point.";
}

TEST(Overapproximation, ConstraintReduction)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/overapproximation/rand_conzono.json";

    ZonoPtr Z = from_json(filename);

    ConZono* Z_cz_ptr = dynamic_cast<ConZono*>(Z.get());
    ASSERT_NE(Z_cz_ptr, nullptr) << "Expected set to be a constrained zonotope for constraint reduction test.";
    Z.release();
    std::unique_ptr<ConZono> Z_conzono (Z_cz_ptr);

    Eigen::Vector<zono_float, -1> p (Z_conzono->get_n());
    p.setZero();
    const auto p_proj = Z_conzono->project_point(p);

    const auto Z_cr = Z_conzono->constraint_reduction();

    EXPECT_TRUE(Z_cr->contains_point(p_proj)) << "Reduced constrained zonotope does not contain projected point.";
}

TEST(Overapproximation, ToZonoApprox)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/overapproximation/rand_conzono.json";

    ZonoPtr Z = from_json(filename);

    ConZono* Z_conzono_ptr = dynamic_cast<ConZono*>(Z.get());
    ASSERT_NE(Z_conzono_ptr, nullptr) << "Expected set to be a constrained zonotope for zonotope approximation test.";
    Z.release();
    std::unique_ptr<ConZono> Z_conzono (Z_conzono_ptr);

    Eigen::Vector<zono_float, -1> p (Z_conzono->get_n());
    p.setZero();
    const auto p_proj = Z_conzono->project_point(p);

    const auto Z_approx = Z_conzono->to_zono_approx();

    EXPECT_TRUE(Z_approx->contains_point(p_proj)) << "Zonotope approximation does not contain projected point.";
}

TEST(Overapproximation, ConvexRelaxation)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set";
    const std::string filename = g_test_data_dir + "/overapproximation/rand_hybzono.json";

    ZonoPtr Z = from_json(filename);

    Eigen::Vector<zono_float, -1> p (Z->get_n());
    p.setZero();
    const auto p_proj = Z->project_point(p);

    const auto Z_relax = Z->convex_relaxation();

    EXPECT_TRUE(Z_relax->contains_point(p_proj)) << "Convex relaxation does not contain projected point.";
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc >= 2) g_test_data_dir = argv[1];
    return RUN_ALL_TESTS();
}
