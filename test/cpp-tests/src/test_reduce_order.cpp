#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

class ReduceOrderTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // regular 16-sided zonotope (8 generators), centered away from the origin
        Z = make_regular_zono_2D(1.0, 16, false, Eigen::Vector<zono_float, 2>(1.0, -2.0));

        // unit directions at 16 evenly spaced angles
        for (int k = 0; k < 16; ++k)
        {
            const zono_float th = 2 * pi * k / 16;
            Eigen::Vector<zono_float, -1> d(2);
            d << std::cos(th), std::sin(th);
            directions.push_back(d);
        }
    }

    std::unique_ptr<Zono> Z;
    std::vector<Eigen::Vector<zono_float, -1>> directions;
};

TEST_F(ReduceOrderTest, ReducedSetContainsOriginal)
{
    std::mt19937 rand_gen(0);
    std::uniform_real_distribution<zono_float> unif(-1, 1);
    const Eigen::Matrix<zono_float, -1, -1> G = Z->get_G().toDense();

    for (const int n_o : {2, 3, 5, 7})
    {
        const auto Zr = Z->reduce_order(n_o);
        EXPECT_EQ(Zr->get_nG(), n_o);

        // random points of Z: x = G xi + c with xi uniform in [-1, 1]^nG
        for (int s = 0; s < 50; ++s)
        {
            Eigen::Vector<zono_float, -1> xi(Z->get_nG());
            for (int k = 0; k < xi.size(); ++k) xi(k) = unif(rand_gen);
            const Eigen::Vector<zono_float, -1> x = G * xi + Z->get_c();
            EXPECT_TRUE(Zr->contains_point(x)) << "n_o = " << n_o << ": reduced set must contain sample "
                                               << x.transpose();
        }
    }
}

TEST_F(ReduceOrderTest, NoReductionNeededOnlySorts)
{
    for (const int n_o : {8, 13})
    {
        const auto Zr = Z->reduce_order(n_o);

        EXPECT_EQ(Zr->get_nG(), Z->get_nG()) << "n_o = " << n_o;
        for (const auto& d : directions)
        {
            EXPECT_NEAR(Zr->support(d), Z->support(d), 1e-9) << "n_o = " << n_o << ": set must be unchanged";
        }

        const Eigen::Matrix<zono_float, -1, -1> Gr = Zr->get_G().toDense();
        for (int k = 1; k < Gr.cols(); ++k)
        {
            EXPECT_GE(Gr.col(k - 1).norm(), Gr.col(k).norm() - 1e-12)
                << "n_o = " << n_o << ": generators must be sorted by decreasing norm";
        }
    }
}

TEST_F(ReduceOrderTest, OrderLessThanDimensionThrows)
{
    EXPECT_THROW(Z->reduce_order(1), std::invalid_argument);
}
