#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;
using namespace ZonoOpt::detail;

namespace
{
    // trivial 2-variable QP: min 0.5*||x||^2 s.t. x1 + x2 = 1
    ADMM_data make_simple_admm_data(const bool use_interval_contractor)
    {
        Eigen::SparseMatrix<zono_float> P(2, 2);
        P.setIdentity();

        Eigen::SparseMatrix<zono_float> A(1, 2);
        A.insert(0, 0) = 1.0;
        A.insert(0, 1) = 1.0;
        A.makeCompressed();

        const Eigen::Vector<zono_float, -1> q = Eigen::Vector<zono_float, -1>::Zero(2);
        Eigen::Vector<zono_float, -1> b(1);
        b << 1.0;
        const Eigen::Vector<zono_float, -1> x_l = -10.0 * Eigen::Vector<zono_float, -1>::Ones(2);
        const Eigen::Vector<zono_float, -1> x_u = 10.0 * Eigen::Vector<zono_float, -1>::Ones(2);

        OptSettings settings;
        settings.use_interval_contractor = use_interval_contractor;
        settings.verbose = false;

        ADMM_data data;
        data.set(P, q, A, b, x_l, x_u, 0.0, settings);
        return data;
    }
}

// clone() should share the ADMM_prob core rather than deep-copying it, while x_box stays per-instance
TEST(AdmmDataSharing, CloneSharesProb)
{
    const ADMM_data d = make_simple_admm_data(false);
    const std::shared_ptr<ADMM_data> c1(d.clone());
    const std::shared_ptr<ADMM_data> c2(c1->clone());

    EXPECT_EQ(d.prob.get(), c1->prob.get());
    EXPECT_EQ(c1->prob.get(), c2->prob.get());
    EXPECT_NE(d.x_box.get(), c1->x_box.get());
}

// q is per-instance (ConZono::do_bounding_box relies on mutating it independently of sibling clones)
TEST(AdmmDataSharing, MutatingQDoesNotAffectSiblingClone)
{
    const ADMM_data d = make_simple_admm_data(false);
    const std::shared_ptr<ADMM_data> c1(d.clone());
    const std::shared_ptr<ADMM_data> c2(c1->clone());

    c1->q(0) = 42.0;
    EXPECT_NE(c1->q(0), c2->q(0));
}

// A_rm is lazily built only when the interval contractor is actually used
TEST(AdmmDataSharing, ARmBuiltOnlyWhenContractorUsed)
{
    {
        const auto data = std::make_shared<ADMM_data>(make_simple_admm_data(false));
        ADMM_solver solver(data);
        solver.solve();
        EXPECT_FALSE(data->prob->A_rm_ready());
    }
    {
        const auto data = std::make_shared<ADMM_data>(make_simple_admm_data(true));
        ADMM_solver solver(data);
        solver.solve();
        EXPECT_TRUE(data->prob->A_rm_ready());
    }
}

// the shared factorization is only valid for the rho it was built with
TEST(AdmmDataSharing, RhoMismatchThrows)
{
    const auto data = std::make_shared<ADMM_data>(make_simple_admm_data(false));
    data->settings.rho *= 2;
    ADMM_solver solver(data);
    EXPECT_THROW(solver.solve(), std::runtime_error);
}
