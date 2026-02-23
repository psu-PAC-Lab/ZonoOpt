#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>

using namespace ZonoOpt;

int main()
{
    // make random conzonos
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

    // take union
    const auto U = union_of_many(CZs);

    // minkowski sum
    const auto Z = minkowski_sum(*U, *U);

    // get number of leaves
    const auto leaves = Z->get_leaves();

    // check number of leaves is correct
    test_assert(leaves.size() == n_CZs * n_CZs, "Expected " + std::to_string(n_CZs * n_CZs) + " leaves, got " + std::to_string(leaves.size()));

    return 0;
}