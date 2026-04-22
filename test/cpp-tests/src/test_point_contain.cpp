#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

static std::string g_test_data_dir;

TEST(PointContain, ContainsInternalPoint)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/point_contain/";

    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "c.txt");
    const Eigen::SparseMatrix<zono_float> A = load_sparse_matrix(test_folder + "A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "b.txt");

    ConZono Z (G, c, A, b);

    const Eigen::Vector<zono_float, -1> x_c = load_vector(test_folder + "x_c.txt");

    EXPECT_TRUE(Z.contains_point(x_c)) << "Expected Z to contain x_c";
}

TEST(PointContain, RejectsExternalPoint)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/point_contain/";

    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "c.txt");
    const Eigen::SparseMatrix<zono_float> A = load_sparse_matrix(test_folder + "A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "b.txt");

    ConZono Z (G, c, A, b);

    const Eigen::Vector<zono_float, -1> x_n = load_vector(test_folder + "x_n.txt");

    EXPECT_FALSE(Z.contains_point(x_n)) << "Expected Z to not contain x_n";
}

TEST(PointContain, SetNotFullDimensional)
{
    const std::string filename = g_test_data_dir + "/point_contain/Z_not_full_dim.json";
    auto Z = from_json(filename);
    auto Z_cr = Z->convex_relaxation();
    auto Z_rr = Z_cr->remove_redundancy();

    // allow QR factorization for rank-deficient A in ADMM
    OptSettings settings;
    settings.rank_deficient_qr_admm = true;

    // project c vector onto set
    const Eigen::Vector<zono_float, -1> c_proj = Z_rr->project_point(Z_rr->get_c(), settings);

    // check that projected point is contained
    EXPECT_TRUE(Z_rr->contains_point(c_proj, settings)) << "Expected convex relaxation to contain projected point";
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc >= 2) g_test_data_dir = argv[1];
    return RUN_ALL_TESTS();
}
