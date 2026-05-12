#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

TEST(ZonoHull, PointsContained)
{
    // create two random zonotopes
    const int n = 10;
    const double density = 1.;
    double val_min = -1.;
    double val_max = 1.;
    std::mt19937 rand_gen(0);
    auto Z1 = random_zono(n, 30, density, val_min, val_max, rand_gen);
    auto Z2 = random_zono(n, 40, density, val_min, val_max, rand_gen);
    
    // random offset
    Eigen::Vector<zono_float, -1> offset(n);
    val_min = -100.;
    val_max = 100.;
    for (int i = 0; i < n; ++i)
    {
        offset(i) = (val_max - val_min) * rand_gen() / rand_gen.max() + val_min;
    }
    Z2 += offset;

    // compute hull
    std::vector<std::shared_ptr<HybZono>> Zs;
    Zs.push_back(std::make_shared<Zono>(Z1));
    Zs.push_back(std::make_shared<Zono>(Z2));

    auto U = convex_hull(Zs, false);

    // check that points along line between centers are in set
    EXPECT_TRUE(U->contains_point(Z1.get_center()));
    EXPECT_TRUE(U->contains_point(Z2.get_center()));
    EXPECT_TRUE(U->contains_point((Z1.get_center() + Z2.get_center()) / 2));
}

TEST(ZonoHull, SeveralSets)
{
    // create several random zonotopes
    const int n = 10;
    const double density = 1.;
    double val_min = -1.;
    double val_max = 1.;
    std::mt19937 rand_gen(0);
    std::vector<std::shared_ptr<HybZono>> Zs;
    for (int i = 0; i < 5; ++i)
    {
        auto Z = random_zono(n, 30, density, val_min, val_max, rand_gen);
        Zs.push_back(std::make_shared<Zono>(Z));
    }

    auto U = convex_hull(Zs, false);

    // check that points in each zonotope are in the hull
    for (const auto& Z : Zs)
    {
        Zono* Z_zono = dynamic_cast<Zono*>(Z.get());
        ASSERT_TRUE(Z_zono != nullptr);
        EXPECT_TRUE(U->contains_point(Z_zono->get_center()));
    }
}

TEST(ZonoHull, InconsistentDimensions)
{
    // create two zonotopes with different dimensions
    const int n1 = 10;
    const int n2 = 15;
    const double density = 1.;
    double val_min = -1.;
    double val_max = 1.;
    std::mt19937 rand_gen(0);
    auto Z1 = random_zono(n1, 30, density, val_min, val_max, rand_gen);
    auto Z2 = random_zono(n2, 40, density, val_min, val_max, rand_gen);

    // put them in a vector
    std::vector<std::shared_ptr<HybZono>> Zs;
    Zs.push_back(std::make_shared<Zono>(Z1));
    Zs.push_back(std::make_shared<Zono>(Z2));

    // try to compute the convex hull - this should fail
    EXPECT_THROW(convex_hull(Zs, false), std::invalid_argument);
}