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

    if (detail::gurobi_available())
    {
        EXPECT_FALSE(Zf.is_empty(GurobiSettings())) << "Expected Zf to be non-empty using Gurobi";
    }
    if (detail::scip_available())
    {
        EXPECT_FALSE(Zf.is_empty(SCIPSettings())) << "Expected Zf to be non-empty using SCIP";
    }
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

    if (detail::gurobi_available())
    {
        EXPECT_TRUE(Zi.is_empty(GurobiSettings())) << "Expected Zi to be empty using Gurobi";
    }
    if (detail::scip_available())
    {
        EXPECT_TRUE(Zi.is_empty(SCIPSettings())) << "Expected Zi to be empty using SCIP";
    }
}

class EmptyConstraintRowTest : public ::testing::Test
{
protected:
    // unit box with constraint rows that have no nonzero coefficients
    // (external solvers only: ADMM requires A to have full row rank)
    Eigen::SparseMatrix<zono_float> G = Eigen::SparseMatrix<zono_float>(2, 2);
    Eigen::SparseMatrix<zono_float> A_one_empty_row = Eigen::SparseMatrix<zono_float>(2, 2); // row 0 empty
    Eigen::SparseMatrix<zono_float> A_all_empty = Eigen::SparseMatrix<zono_float>(1, 2);
    Eigen::Vector<zono_float, -1> c = Eigen::Vector<zono_float, -1>::Zero(2);

    void SetUp() override
    {
        G.setIdentity();
        A_one_empty_row.insert(1, 0) = 1;
    }

    static void expect_is_empty(const ConZono& Z, const bool expected, const std::string& name)
    {
        if (detail::gurobi_available())
        {
            EXPECT_EQ(Z.is_empty(GurobiSettings()), expected) << name << " using Gurobi";
        }
        if (detail::scip_available())
        {
            EXPECT_EQ(Z.is_empty(SCIPSettings()), expected) << name << " using SCIP";
        }
    }
};

TEST_F(EmptyConstraintRowTest, ZeroRightHandSideIsNotEmpty)
{
    // 0 = 0 is always satisfied
    expect_is_empty(ConZono(G, c, A_one_empty_row, Eigen::Vector<zono_float, -1>::Zero(2)), false, "one empty row");
    expect_is_empty(ConZono(G, c, A_all_empty, Eigen::Vector<zono_float, -1>::Zero(1)), false, "all rows empty");
}

TEST_F(EmptyConstraintRowTest, NonzeroRightHandSideIsEmpty)
{
    // 0 = 1 can never be satisfied
    Eigen::Vector<zono_float, -1> b(2);
    b << 1, 0;
    expect_is_empty(ConZono(G, c, A_one_empty_row, b), true, "one empty row");
    expect_is_empty(ConZono(G, c, A_all_empty, Eigen::Vector<zono_float, -1>::Ones(1)), true, "all rows empty");
}
