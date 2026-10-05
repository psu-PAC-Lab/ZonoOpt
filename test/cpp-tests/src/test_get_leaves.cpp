#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>

using namespace ZonoOpt;

TEST(GetLeaves, LeavesCount)
{
    constexpr int n_CZs = 20;
    constexpr int n = 10;
    constexpr int nV = 2 * n;

    srand(0);
    std::vector<std::shared_ptr<HybZono>> CZs;
    for (int i = 0; i < n_CZs; i++)
    {
        Eigen::Matrix<zono_float, -1, -1> V = Eigen::Matrix<zono_float, -1, -1>::Random(nV, n);
        CZs.push_back(vrep_2_conzono(V));
    }

    const auto U = union_of_many(CZs);
    const auto Z = minkowski_sum(*U, *U);
    const auto leaves = Z->get_leaves();

    EXPECT_EQ(leaves.size(), static_cast<size_t>(n_CZs * n_CZs))
        << "Expected " << n_CZs * n_CZs << " leaves, got " << leaves.size();

    if (detail::gurobi_available())
    {
        const auto leaves_grb = Z->get_leaves(GetLeavesParams(), GurobiSettings());
        EXPECT_EQ(leaves_grb.size(), static_cast<size_t>(n_CZs * n_CZs))
            << "Expected " << n_CZs * n_CZs << " leaves using Gurobi, got " << leaves_grb.size();
    }
    if (detail::scip_available())
    {
        const auto leaves_scip = Z->get_leaves(GetLeavesParams(), SCIPSettings());
        EXPECT_EQ(leaves_scip.size(), static_cast<size_t>(n_CZs * n_CZs))
            << "Expected " << n_CZs * n_CZs << " leaves using SCIP, got " << leaves_scip.size();
    }
}

TEST(GetLeaves, ScipSolutionSatisfiesConstraints)
{
    if (!detail::scip_available())
        GTEST_SKIP() << "SCIP not available";

    // union of two offset boxes, so the continuous factors of a solution are nonzero
    std::vector<std::shared_ptr<Zono>> Zs;
    Zs.push_back(make_regular_zono_2D(1, 4, false, Eigen::Vector<zono_float, 2>(3, 0)));
    Zs.push_back(make_regular_zono_2D(1, 4, false, Eigen::Vector<zono_float, 2>(-3, 0)));
    const auto Z = zono_union_2_hybzono(Zs);

    std::shared_ptr<OptSolution> sol;
    const auto leaves = Z->get_leaves(GetLeavesParams(), SCIPSettings(), &sol);
    EXPECT_EQ(leaves.size(), 2u);

    // get_leaves solves min 0.5 xi^T xi s.t. A xi = b
    ASSERT_TRUE(sol);
    EXPECT_LT((Z->get_A() * sol->z - Z->get_b()).cwiseAbs().maxCoeff(), 1e-4)
        << "SCIP multisol solution must satisfy A xi = b";
    EXPECT_NEAR(sol->J, 0.5 * sol->z.squaredNorm(), 1e-4) << "SCIP multisol objective must match its solution";
}
