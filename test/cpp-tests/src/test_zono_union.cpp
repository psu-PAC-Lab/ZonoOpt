#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

namespace
{
    std::shared_ptr<Zono> make_zono(const Eigen::Matrix<zono_float, -1, -1>& G, const Eigen::Vector<zono_float, -1>& c)
    {
        return std::make_shared<Zono>(Eigen::SparseMatrix<zono_float>(G.sparseView()), c);
    }

    Eigen::Vector<zono_float, -1> vec(std::initializer_list<zono_float> vals)
    {
        Eigen::Vector<zono_float, -1> v(static_cast<Eigen::Index>(vals.size()));
        Eigen::Index k = 0;
        for (const zono_float x : vals) v(k++) = x;
        return v;
    }
}

TEST(ZonoUnion, DistinctGenerators)
{
    // unit box at the origin union a diamond centered at (5, 0)
    Eigen::Matrix<zono_float, -1, -1> G1(2, 2), G2(2, 2);
    G1 << 1, 0,
          0, 1;
    G2 << 1, 1,
          1, -1;
    std::vector<std::shared_ptr<Zono>> Zs = {make_zono(G1, vec({0, 0})), make_zono(G2, vec({5, 0}))};
    const auto U = zono_union_2_hybzono(Zs);

    EXPECT_TRUE(U->contains_point(vec({0.5, -0.5})));
    EXPECT_TRUE(U->contains_point(vec({-0.9, 0.9})));
    EXPECT_TRUE(U->contains_point(vec({6.5, 0})));
    EXPECT_TRUE(U->contains_point(vec({5, 1.5})));
    EXPECT_FALSE(U->contains_point(vec({2.5, 0})));
    EXPECT_FALSE(U->contains_point(vec({6.5, 1})));
    EXPECT_FALSE(U->contains_point(vec({0, 1.5})));
}

TEST(ZonoUnion, DuplicatedGeneratorsBetweenZonotopes)
{
    // unit box at the origin union a box [4, 6] x [-2, 2]; both use the generator [1, 0]^T, and the
    // first also uses [0, 1]^T
    Eigen::Matrix<zono_float, -1, -1> G1(2, 2), G2(2, 2);
    G1 << 1, 0,
          0, 1;
    G2 << 1, 0,
          0, 2;
    std::vector<std::shared_ptr<Zono>> Zs = {make_zono(G1, vec({0, 0})), make_zono(G2, vec({5, 0}))};
    const auto U = zono_union_2_hybzono(Zs);

    // 3 unique generators, each with one factor and one slack factor
    EXPECT_EQ(U->get_nGc(), 6);

    EXPECT_TRUE(U->contains_point(vec({0.5, -0.5})));
    EXPECT_TRUE(U->contains_point(vec({-0.9, 0.9})));
    EXPECT_TRUE(U->contains_point(vec({5.5, 1.5})));
    EXPECT_TRUE(U->contains_point(vec({4.5, -1.8})));
    EXPECT_FALSE(U->contains_point(vec({2.5, 0})));
    EXPECT_FALSE(U->contains_point(vec({0, 1.5})));
    EXPECT_FALSE(U->contains_point(vec({6.5, 0})));
}

TEST(ZonoUnion, DuplicatedGeneratorInOneZonotope)
{
    // G = [[1, 1, 0], [0, 0, 1]] is the box [-2, 2] x [-1, 1]; both equal columns must be kept
    Eigen::Matrix<zono_float, -1, -1> G(2, 3);
    G << 1, 1, 0,
         0, 0, 1;
    std::vector<std::shared_ptr<Zono>> Zs = {make_zono(G, vec({0, 0}))};
    const auto U = zono_union_2_hybzono(Zs);

    EXPECT_TRUE(U->contains_point(vec({1.5, 0.5})));
    EXPECT_TRUE(U->contains_point(vec({-1.5, -0.5})));
    EXPECT_FALSE(U->contains_point(vec({2.5, 0})));
}
