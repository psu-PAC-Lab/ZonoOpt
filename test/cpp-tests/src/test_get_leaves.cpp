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
        const auto leaves_grb = Z->get_leaves(false, GurobiSettings());
        EXPECT_EQ(leaves_grb.size(), static_cast<size_t>(n_CZs * n_CZs))
            << "Expected " << n_CZs * n_CZs << " leaves using Gurobi, got " << leaves_grb.size();
    }
    if (detail::scip_available())
    {
        const auto leaves_scip = Z->get_leaves(false, SCIPSettings());
        EXPECT_EQ(leaves_scip.size(), static_cast<size_t>(n_CZs * n_CZs))
            << "Expected " << n_CZs * n_CZs << " leaves using SCIP, got " << leaves_scip.size();
    }
}