#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

namespace
{
    Eigen::SparseMatrix<zono_float> from_triplets(const int rows, const int cols,
                                                  const std::vector<Eigen::Triplet<zono_float>>& triplets)
    {
        Eigen::SparseMatrix<zono_float> M(rows, cols);
        M.setFromTriplets(triplets.begin(), triplets.end());
        return M;
    }

    // cost 0.5 x^T x, i.e. P = I, q = 0
    void expect_nan_point(const HybZono& Z, const std::string& name)
    {
        Eigen::SparseMatrix<zono_float> P(Z.get_n(), Z.get_n());
        P.setIdentity();
        const Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(Z.get_n());

        std::shared_ptr<OptSolution> sol;
        const Eigen::Vector<zono_float, -1> x = Z.optimize_over(P, q, 0, get_default_solver_settings(), &sol);

        ASSERT_TRUE(sol && sol->infeasible) << name << ": expected the problem to be reported infeasible";
        EXPECT_EQ(x.size(), Z.get_n()) << name << ": result must have the set dimension n, not nG";
        EXPECT_TRUE(x.array().isNaN().all()) << name << ": infeasible result must be all NaN";
    }
}

TEST(OptimizeOver, InfeasibleConZonoReturnsNaN)
{
    // n = 2, nG = 3; the constraint xi_0 = 5 cannot hold for xi in [-1, 1]^3
    const auto G = from_triplets(2, 3, {{0, 0, 1}, {1, 1, 1}, {0, 2, 1}, {1, 2, 1}});
    const auto A = from_triplets(1, 3, {{0, 0, 1}});
    const ConZono Z(G, Eigen::Vector<zono_float, -1>::Zero(2), A, Eigen::Vector<zono_float, -1>::Constant(1, 5));

    expect_nan_point(Z, "ConZono");
}

TEST(OptimizeOver, InfeasibleHybZonoReturnsNaN)
{
    // n = 2, nGc = 2, nGb = 1; the constraint xi_c0 = 5 cannot hold for xi_c in [-1, 1]^2
    const auto Gc = from_triplets(2, 2, {{0, 0, 1}, {1, 1, 1}});
    const auto Gb = from_triplets(2, 1, {{0, 0, 1}});
    const auto Ac = from_triplets(1, 2, {{0, 0, 1}});
    const auto Ab = from_triplets(1, 1, {});
    const HybZono Z(Gc, Gb, Eigen::Vector<zono_float, -1>::Zero(2), Ac, Ab,
                    Eigen::Vector<zono_float, -1>::Constant(1, 5));

    expect_nan_point(Z, "HybZono");
}

TEST(OptimizeOver, TimeoutIsNotInfeasible)
{
    // feasible hybrid zonotope: Cartesian product of unions of offset regular zonotopes
    std::vector<std::shared_ptr<Zono>> Zs;
    for (int i = 0; i < 20; ++i)
    {
        Eigen::Vector<zono_float, 2> c;
        c << static_cast<zono_float>(3 * i), static_cast<zono_float>(i % 2);
        Zs.push_back(make_regular_zono_2D(1, 8, false, c));
    }
    std::unique_ptr<HybZono> U = zono_union_2_hybzono(Zs);
    std::unique_ptr<HybZono> Z = cartesian_product(*U, *U);
    for (int i = 0; i < 3; ++i)
        Z = cartesian_product(*Z, *U);

    // time limit too short to find a solution
    OptSettings settings;
    settings.t_max = 1e-9;

    Eigen::SparseMatrix<zono_float> P(Z->get_n(), Z->get_n());
    P.setIdentity();
    const Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(Z->get_n());

    std::shared_ptr<OptSolution> sol;
    const Eigen::Vector<zono_float, -1> x = Z->optimize_over(P, q, 0, settings, &sol);

    ASSERT_TRUE(sol);
    EXPECT_FALSE(sol->infeasible) << "a timeout must not be reported as proven infeasibility";
    EXPECT_FALSE(sol->converged);
    EXPECT_EQ(x.size(), Z->get_n());
}

class ScipQpTest : public ::testing::Test
{
protected:
    Eigen::SparseMatrix<zono_float> P;
    Eigen::Vector<zono_float, -1> q;

    void SetUp() override
    {
        if (!detail::scip_available())
            GTEST_SKIP() << "SCIP not available";
    }

    // returns {J_admm, J_scip}, the objectives at the ADMM and SCIP solutions
    std::pair<zono_float, zono_float> solve_both(const HybZono& Z, std::shared_ptr<OptSolution>* scip_sol) const
    {
        OptSettings admm_settings;
        admm_settings.eps_prim = static_cast<zono_float>(1e-6);
        admm_settings.eps_dual = static_cast<zono_float>(1e-6);
        const Eigen::Vector<zono_float, -1> x_admm = Z.optimize_over(P, q, 0, admm_settings);
        const zono_float J_admm = 0.5 * x_admm.dot(P * x_admm) + q.dot(x_admm);

        const Eigen::Vector<zono_float, -1> x_scip = Z.optimize_over(P, q, 0, SCIPSettings(), scip_sol);
        const zono_float J_scip = 0.5 * x_scip.dot(P * x_scip) + q.dot(x_scip);
        return {J_admm, J_scip};
    }
};

TEST_F(ScipQpTest, OriginInfeasible)
{
    // random feasible constrained zonotopes whose constraints exclude xi = 0
    std::mt19937 gen(0);
    for (int trial = 0; trial < 5; ++trial)
    {
        constexpr int n = 4, nG = 12, nC = 4;
        const auto G = random_sparse_matrix(n, nG, 0.6, -1, 1, gen);
        const auto A = random_sparse_matrix(nC, nG, 0.6, -1, 1, gen);
        const Eigen::Vector<zono_float, -1> b = A * random_vector(nG, -0.5, 0.5, gen);
        const ConZono Z(G, Eigen::Vector<zono_float, -1>::Zero(n), A, b);

        P.resize(n, n);
        P.setIdentity();
        q = random_vector(n, -5, 5, gen);
        std::shared_ptr<OptSolution> sol;
        const auto [J_admm, J_scip] = solve_both(Z, &sol);
        EXPECT_NEAR(J_scip, J_admm, 1e-3 * (1 + std::abs(J_admm))) << "random ConZono #" << trial;
    }
}

TEST_F(ScipQpTest, OriginFeasible)
{
    // zonotope with no constraints, so xi = 0 is feasible but not optimal
    P.resize(2, 2);
    P.setIdentity();
    q.resize(2);
    q << 3, -2;
    std::shared_ptr<OptSolution> sol;
    const auto [J_admm, J_scip] = solve_both(*make_regular_zono_2D(1, 8), &sol);
    ASSERT_TRUE(sol);
    if (sol->converged)
    {
        EXPECT_NEAR(J_scip, J_admm, 1e-3 * (1 + std::abs(J_admm)))
            << "SCIP reported convergence to a non-optimal point";
    }
}
