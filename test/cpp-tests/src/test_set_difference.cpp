#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

namespace
{
    // box domain D = [x_lb, x_ub] x [y_lb, y_ub] as a zonotope
    std::unique_ptr<Zono> make_domain(const zono_float x_lb, const zono_float x_ub,
                                      const zono_float y_lb, const zono_float y_ub)
    {
        Eigen::Vector<zono_float, -1> lb(2), ub(2);
        lb << x_lb, y_lb;
        ub << x_ub, y_ub;
        return interval_2_zono(Box(lb, ub));
    }

    // Checks D \ Z on a grid over D for a 2D zonotope Z = G xi with invertible G,
    // classifying points by the exact factor norm r = ||G^-1 p||_inf.
    void check_set_diff_2d(const Eigen::Matrix<zono_float, 2, 2>& Gd)
    {
        const auto D = make_domain(-3, 3, -1, 1);
        Zono Z(Gd.sparseView(), Eigen::Vector<zono_float, -1>::Zero(2));
        const auto D_minus_Z = set_diff(*D, Z, 10);

        const auto lu = Gd.fullPivLu();
        int n_out = 0, n_in = 0;
        for (zono_float x = -2.75; x <= 2.751; x += 0.5)
        {
            for (zono_float y = -0.75; y <= 0.751; y += 0.25)
            {
                const Eigen::Vector<zono_float, 2> p(x, y);
                const zono_float r = lu.solve(p).cwiseAbs().maxCoeff();
                if (r > 1.05)
                {
                    ++n_out;
                    EXPECT_TRUE(D_minus_Z->contains_point(p))
                        << "point (" << x << ", " << y << ") is in D and outside Z, so it must be in D \\ Z";
                }
                else if (r < 0.95)
                {
                    ++n_in;
                    EXPECT_FALSE(D_minus_Z->contains_point(p))
                        << "point (" << x << ", " << y << ") is in the interior of Z, so it must not be in D \\ Z";
                }
            }
        }
        ASSERT_GT(n_out, 0);
        ASSERT_GT(n_in, 0);

        // outside the domain
        EXPECT_FALSE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(5, 0)));
        EXPECT_FALSE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(0, 3)));
    }
}

TEST(SetDifference, CancellingGenerators)
{
    Eigen::Matrix<zono_float, 2, 2> G;
    G << -4, 1,
          1, 0;
    check_set_diff_2d(G);
}

TEST(SetDifference, PositiveDiagonalGenerators)
{
    Eigen::Matrix<zono_float, 2, 2> G;
    G << 2, 0,
         0, 0.5;
    check_set_diff_2d(G);
}

TEST(SetDifference, RedundantConstraints)
{
    // unit box written with an extra generator that is fixed to zero by two linearly dependent constraints
    const std::vector<Eigen::Triplet<zono_float>> trip_G = {{0, 0, 1}, {1, 1, 1}, {0, 2, 1}};
    Eigen::SparseMatrix<zono_float> G(2, 3);
    G.setFromTriplets(trip_G.begin(), trip_G.end());
    const std::vector<Eigen::Triplet<zono_float>> trip_A = {{0, 2, 1}, {1, 2, 2}};
    Eigen::SparseMatrix<zono_float> A(2, 3);
    A.setFromTriplets(trip_A.begin(), trip_A.end());
    ConZono Z(G, Eigen::Vector<zono_float, -1>::Zero(2), A, Eigen::Vector<zono_float, -1>::Zero(2));

    const auto D = make_domain(-3, 3, -3, 3);
    const auto D_minus_Z = set_diff(*D, Z, 10);

    EXPECT_TRUE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(1.5, 0)));
    EXPECT_TRUE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(-2.5, 2)));
    EXPECT_FALSE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(0, 0)));
    EXPECT_FALSE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(0.5, -0.5)));
    EXPECT_FALSE(D_minus_Z->contains_point(Eigen::Vector<zono_float, 2>(4, 0)));
}

TEST(SetDifference, NotFullDimensionalThrows)
{
    // a segment in 2D: [G; A] = [1; 0] does not have full row rank
    const std::vector<Eigen::Triplet<zono_float>> trip_G = {{0, 0, 1}};
    Eigen::SparseMatrix<zono_float> G(2, 1);
    G.setFromTriplets(trip_G.begin(), trip_G.end());
    Zono Z(G, Eigen::Vector<zono_float, -1>::Zero(2));

    const auto D = make_domain(-3, 3, -3, 3);
    EXPECT_THROW(set_diff(*D, Z, 10), std::invalid_argument);
}
