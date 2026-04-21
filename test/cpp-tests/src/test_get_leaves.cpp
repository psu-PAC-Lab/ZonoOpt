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
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
