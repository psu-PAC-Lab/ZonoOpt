#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

const std::string g_test_data_dir = TEST_DATA_DIR;

TEST(IsEmpty, FeasibleIsNotEmpty)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/is_empty/";

    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "f_G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "f_c.txt");
    const Eigen::SparseMatrix<zono_float> A = load_sparse_matrix(test_folder + "f_A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "f_b.txt");

    const ConZono Zf (G, c, A, b);

    EXPECT_FALSE(Zf.is_empty()) << "Expected Zf to be non-empty";
}

TEST(IsEmpty, InfeasibleIsEmpty)
{
    ASSERT_FALSE(g_test_data_dir.empty()) << "Test data directory not set (pass path as argv[1])";
    const std::string test_folder = g_test_data_dir + "/is_empty/";

    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "i_G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "i_c.txt");
    const Eigen::SparseMatrix<zono_float> A = load_sparse_matrix(test_folder + "i_A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "i_b.txt");

    const ConZono Zi (G, c, A, b);

    EXPECT_TRUE(Zi.is_empty()) << "Expected Zi to be empty";
}