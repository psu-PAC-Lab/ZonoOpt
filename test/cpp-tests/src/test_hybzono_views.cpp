#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

namespace
{
    using Dense = Eigen::Matrix<zono_float, -1, -1>;

    Dense dense(const Eigen::SparseMatrix<zono_float>& M)
    {
        return Dense(M);
    }

    // check that Gc/Gb/Ac/Ab getters are consistent with G = [Gc, Gb], A = [Ac, Ab]
    void expect_consistent(const HybZono& Z)
    {
        const Dense Gc = dense(Z.get_Gc()), Gb = dense(Z.get_Gb()), G = dense(Z.get_G());
        const Dense Ac = dense(Z.get_Ac()), Ab = dense(Z.get_Ab()), A = dense(Z.get_A());

        ASSERT_EQ(Gc.cols(), Z.get_nGc());
        ASSERT_EQ(Gb.cols(), Z.get_nGb());
        ASSERT_EQ(G.cols(), Z.get_nG());
        ASSERT_EQ(Ac.cols(), Z.get_nGc());
        ASSERT_EQ(Ab.cols(), Z.get_nGb());
        ASSERT_EQ(A.rows(), Z.get_nC());
        ASSERT_EQ(Ac.rows(), Z.get_nC());
        ASSERT_EQ(Ab.rows(), Z.get_nC());

        Dense G_cat(G.rows(), G.cols()), A_cat(A.rows(), A.cols());
        G_cat << Gc, Gb;
        A_cat << Ac, Ab;
        EXPECT_EQ(G_cat, G);
        EXPECT_EQ(A_cat, A);
    }

    HybZono make_hybzono(const int n, const int nGc, const int nGb, const int nC, std::mt19937& gen)
    {
        auto rnd = [&](const int r, const int c)
        {
            Dense M(r, c);
            for (int i = 0; i < r; ++i)
                for (int j = 0; j < c; ++j)
                    M(i, j) = static_cast<zono_float>(static_cast<int>(gen() % 5) - 2);
            return Eigen::SparseMatrix<zono_float>(M.sparseView());
        };
        Eigen::Vector<zono_float, -1> c = Eigen::Vector<zono_float, -1>::Random(n);
        Eigen::Vector<zono_float, -1> b = Eigen::Vector<zono_float, -1>::Zero(nC);
        return HybZono(rnd(n, nGc), rnd(n, nGb), c, rnd(nC, nGc), rnd(nC, nGb), b, false);
    }
}

TEST(HybZonoViews, GettersRoundTrip)
{
    std::mt19937 gen(0);
    const Dense Gc = Dense::Random(3, 4), Gb = Dense::Random(3, 2);
    const Dense Ac = Dense::Random(2, 4), Ab = Dense::Random(2, 2);
    const Eigen::Vector<zono_float, -1> c = Eigen::Vector<zono_float, -1>::Random(3);
    const Eigen::Vector<zono_float, -1> b = Eigen::Vector<zono_float, -1>::Random(2);

    HybZono Z(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);
    expect_consistent(Z);
    EXPECT_EQ(dense(Z.get_Gc()), Gc);
    EXPECT_EQ(dense(Z.get_Gb()), Gb);
    EXPECT_EQ(dense(Z.get_Ac()), Ac);
    EXPECT_EQ(dense(Z.get_Ab()), Ab);
}

TEST(HybZonoViews, EdgeCases)
{
    std::mt19937 gen(1);
    for (const auto& [nGc, nGb, nC] : std::vector<std::tuple<int, int, int>>{
             {0, 3, 2}, {3, 0, 2}, {3, 2, 0}, {0, 3, 0}, {0, 0, 0}})
    {
        const HybZono Z = make_hybzono(3, nGc, nGb, nC, gen);
        expect_consistent(Z);
    }
}

TEST(HybZonoViews, ConvertFormRoundTrip)
{
    std::mt19937 gen(2);
    HybZono Z = make_hybzono(3, 4, 3, 2, gen);
    const Dense G0 = dense(Z.get_G()), A0 = dense(Z.get_A());
    const Eigen::Vector<zono_float, -1> c0 = Z.get_c(), b0 = Z.get_b();

    Z.convert_form();
    expect_consistent(Z);
    EXPECT_TRUE(Z.is_0_1_form());
    EXPECT_LT((dense(Z.get_G()) - 2 * G0).norm(), 1e-12);
    EXPECT_LT((dense(Z.get_A()) - 2 * A0).norm(), 1e-12);

    Z.convert_form();
    expect_consistent(Z);
    EXPECT_FALSE(Z.is_0_1_form());
    EXPECT_LT((dense(Z.get_G()) - G0).norm(), 1e-12);
    EXPECT_LT((dense(Z.get_A()) - A0).norm(), 1e-12);
    EXPECT_LT((Z.get_c() - c0).norm(), 1e-12);
    EXPECT_LT((Z.get_b() - b0).norm(), 1e-12);
}

TEST(HybZonoViews, ConsistentAfterRemoveRedundancy)
{
    // generator 1 (continuous) and generator 0 (binary) are unused
    Dense Gc(2, 3), Gb(2, 2), Ac(1, 3), Ab(1, 2);
    Gc << 1, 0, 2, 0, 0, 1;
    Gb << 0, 1, 0, 1;
    Ac << 1, 0, 1;
    Ab << 0, 1;
    Eigen::Vector<zono_float, -1> c = Eigen::Vector<zono_float, -1>::Zero(2);
    Eigen::Vector<zono_float, -1> b(1);
    b << 0.5;

    HybZono Z(Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);
    expect_consistent(Z);

    const auto Zr = Z.remove_redundancy();
    ASSERT_NE(Zr, nullptr);
    expect_consistent(*Zr);
}

TEST(HybZonoViews, ConsistentAfterAffineMapAndProduct)
{
    std::mt19937 gen(3);
    const HybZono Z1 = make_hybzono(3, 3, 2, 2, gen);
    HybZono Z2 = make_hybzono(2, 2, 1, 1, gen);

    const Eigen::SparseMatrix<zono_float> R = Dense::Random(4, 3).sparseView();
    const auto Zm = R * Z1;
    expect_consistent(*Zm);
    EXPECT_LT((dense(Zm->get_Gc()) - dense(R * Z1.get_Gc())).norm(), 1e-12);
    EXPECT_LT((dense(Zm->get_Gb()) - dense(R * Z1.get_Gb())).norm(), 1e-12);

    const auto Zp = Z1 * Z2;
    expect_consistent(*Zp);
    EXPECT_EQ(Zp->get_nGc(), Z1.get_nGc() + Z2.get_nGc());
    EXPECT_EQ(Zp->get_nGb(), Z1.get_nGb() + Z2.get_nGb());
}

TEST(HybZonoViews, ConsistentAfterIntersectionAndConstrain)
{
    std::mt19937 gen(4);
    const HybZono Z1 = make_hybzono(3, 3, 2, 2, gen);
    HybZono Z2 = make_hybzono(3, 2, 1, 1, gen);
    HybZono Z1m = Z1;

    // intersection: continuous and binary factors from both sets, no new factors
    const auto Zi = Z1m & Z2;
    expect_consistent(*Zi);
    EXPECT_EQ(Zi->get_nGc(), Z1.get_nGc() + Z2.get_nGc());
    EXPECT_EQ(Zi->get_nGb(), Z1.get_nGb() + Z2.get_nGb());
    EXPECT_EQ(Zi->get_nC(), Z1.get_nC() + Z2.get_nC() + Z1.get_n());
    EXPECT_EQ(dense(Zi->get_Gc()).leftCols(Z1.get_nGc()), dense(Z1.get_Gc()));
    EXPECT_EQ(dense(Zi->get_Gb()).leftCols(Z1.get_nGb()), dense(Z1.get_Gb()));

    // halfspace intersection adds one slack continuous factor per constraint
    Eigen::SparseMatrix<zono_float> H = Dense::Random(2, 3).sparseView();
    Eigen::Vector<zono_float, -1> f = Eigen::Vector<zono_float, -1>::Ones(2);
    const auto Zh = halfspace_intersection(Z1m, H, f);
    expect_consistent(*Zh);
    EXPECT_EQ(Zh->get_nGc(), Z1.get_nGc() + 2);
    EXPECT_EQ(Zh->get_nGb(), Z1.get_nGb());
    EXPECT_EQ(Zh->get_nC(), Z1.get_nC() + 2);

    // equality constraint adds no slack
    const auto Ze = constrain(Z1m, H, f, '=');
    expect_consistent(*Ze);
    EXPECT_EQ(Ze->get_nGc(), Z1.get_nGc());
    EXPECT_EQ(Ze->get_nC(), Z1.get_nC() + 2);
}
