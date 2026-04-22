#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

const std::string g_test_data_dir = TEST_DATA_DIR;

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
    // non-full dimensional constrained zonotope
    Eigen::Matrix<zono_float, 2, 3> G;
    G << 1., 0., 1.,
         0., 1., 1.;
    Eigen::Vector<zono_float, 2> c;
    c << 1., 2.;
    Eigen::Matrix<zono_float, 1, 3> A;
    A << 1., 0., 1.;
    Eigen::Vector<zono_float, 1> b;
    b << 0.5;
    ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // point inside set
    Eigen::Vector<zono_float, -1> x_inside (2);
    x_inside << 1.5, 2.;

    // point outside set
    Eigen::Vector<zono_float, -1> x_outside (2);
    x_outside << 1.7, 2.;

    // allow QR factorization for rank-deficient A in ADMM
    OptSettings settings;
    settings.rank_deficient_qr_admm = true;

    // check that projected point is contained
    EXPECT_TRUE(Z.contains_point(x_inside, settings)) << "Expected convex relaxation to contain projected point";
    EXPECT_FALSE(Z.contains_point(x_outside, settings)) << "Expected convex relaxation to not contain projected point";
}