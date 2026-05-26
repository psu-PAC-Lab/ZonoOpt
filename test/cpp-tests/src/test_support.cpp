#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

const std::string g_test_data_dir = TEST_DATA_DIR;

TEST(Support, ConZonoSupportValue)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/support/";

    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "c.txt");
    const Eigen::SparseMatrix<zono_float> A = load_sparse_matrix(test_folder + "A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "b.txt");

    ConZono Z (G, c, A, b);

    Eigen::Vector<zono_float, -1> d = load_vector(test_folder + "d.txt");
    const zono_float s_expected = load_vector(test_folder + "sup.txt")(0);

    const zono_float s = Z.support(d);

    EXPECT_NEAR(s / s_expected, 1.0, 5e-2) << "Support value does not match expected value";

    if (detail::gurobi_available())
    {
        const zono_float s_grb = Z.support(d, GurobiSettings());
        EXPECT_NEAR(s_grb / s_expected, 1.0, 5e-2) << "Support value using Gurobi does not match expected value";
    }
    if (detail::scip_available())
    {
        const zono_float s_scip = Z.support(d, SCIPSettings());
        EXPECT_NEAR(s_scip / s_expected, 1.0, 5e-2) << "Support value using SCIP does not match expected value";
    }
}

TEST(Support, ZonoConZonoConsistency)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/support/";

    Eigen::Vector<zono_float, -1> d(5);
    d << 1., -1., 0.5, -0.5, 1.5;

    OptSettings settings;
    settings.eps_prim = 1e-3;
    settings.eps_dual = 1e-3;
    settings.rho = 1.;

    for (unsigned int i=0; i<10; ++i)
    {
        const std::string filename = test_folder + "zono_" + std::to_string(i) + ".json";
        auto Zjson = from_json(filename);
        Zono* Zz_tmp = dynamic_cast<Zono*>(Zjson.get());
        ASSERT_NE(Zz_tmp, nullptr) << "Loaded set from " << filename << " is not a Zono.";

        Zjson.release();
        std::unique_ptr<Zono> Zz (Zz_tmp);

        ConZono Zc(Zz->get_G(), Zz->get_c(), Zz->get_A(), Zz->get_b());

        auto sol1 = std::make_shared<OptSolution>();
        auto sol2 = std::make_shared<OptSolution>();
        const zono_float s1 = Zz->support(d, settings, &sol1);
        const zono_float s2 = Zc.support(d, settings, &sol2);

        std::stringstream ss;
        ss << "Zono and ConZono support values do not match"
           << "\nZono support: " << s1 << "\nConZono support: " << s2
           << "\nZ: " << *Zz << ", Zc: " << Zc
           << "\nsol1: " << sol1->print() << "\nsol2: " << sol2->print();
        EXPECT_NEAR(s1, s2, 1e-2) << ss.str();

        ss.str("");
        const zono_float applied1 = d.dot(Zz->get_G() * sol1->z + Zz->get_c());
        const zono_float applied2 = d.dot(Zc.get_G() * sol2->z + Zc.get_c());
        ss << "Zono and ConZono solutions do not produce same support value when applied to generators and center"
           << "\nZono support: " << applied1 << "\nConZono support: " << applied2;
        EXPECT_NEAR(applied1, applied2, 1e-2) << ss.str();

        if (detail::gurobi_available())
        {
            const zono_float s1g = Zz->support(d, GurobiSettings());
            const zono_float s2g = Zc.support(d, GurobiSettings());
            EXPECT_NEAR(s1g, s2g, 1e-2)
                << "Zono and ConZono support values disagree using Gurobi: " << s1g << " vs " << s2g;
        }
        if (detail::scip_available())
        {
            const zono_float s1s = Zz->support(d, SCIPSettings());
            const zono_float s2s = Zc.support(d, SCIPSettings());
            EXPECT_NEAR(s1s, s2s, 1e-2)
                << "Zono and ConZono support values disagree using SCIP: " << s1s << " vs " << s2s;
        }
    }
}