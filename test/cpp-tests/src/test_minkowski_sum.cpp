#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

int main(int argc, char* argv[])
{
    // input: directory where unit test data resides
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <test data folder>" << std::endl;
        return 1;
    }
    std::string test_folder = argv[1];
    test_folder += "/minkowski_sum/";

    // build hybzono from vrep
    std::vector<Eigen::Matrix<zono_float, -1, -1>> V_polys;
    Eigen::Matrix<zono_float, 4, 2> V1;
    V1 << 5.566, 5.896,
        4.044, 5.498,
        5.32, 3.909,
        5.599, 4.082;
    V_polys.push_back(V1);
    Eigen::Matrix<zono_float, 3, 2> V2;
    V2 << 0.049, 6.05,
        -0.248, 3.881,
        0.617, 3.981;
    V_polys.push_back(V2);
    Eigen::Matrix<zono_float, 3, 2> V3;
    V3 << 5.481, 0.911,
        4.937, 1.183,
        5.199, -1.001;
    V_polys.push_back(V3);
    Eigen::Matrix<zono_float, 4, 2> V4;
    V4 << 3.447, 3.207,
        2.853, 3.552,
        3.341, 1.914,
        3.656, 2.397;
    V_polys.push_back(V4);

    const auto Z1 = vrep_2_hybzono(V_polys);

    // zonotope
    Eigen::Matrix<zono_float, 2, 3> G;
    G << 0.5*std::sqrt(3), 0.5, 0.5*std::sqrt(3),
         0.25, 0, -0.25;
    Eigen::Vector2d c(-2.0, 1.0);
    Zono Z2 (G.sparseView(), c);

    // minkowski sum
    const auto Z = minkowski_sum(*Z1, Z2);
    if (Z->is_0_1_form())
        Z->convert_form();

    // expected result
    const Eigen::SparseMatrix<zono_float> Gc_expected = load_sparse_matrix(test_folder + "Gc.txt");
    const Eigen::SparseMatrix<zono_float> Gb_expected = load_sparse_matrix(test_folder + "Gb.txt");
    const Eigen::Vector<zono_float, -1> c_expected = load_vector(test_folder + "c.txt");
    Eigen::SparseMatrix<zono_float> Ac_expected = load_sparse_matrix(test_folder + "Ac.txt");
    Eigen::SparseMatrix<zono_float> Ab_expected = load_sparse_matrix(test_folder + "Ab.txt");
    Eigen::Vector<zono_float, -1> b_expected = load_vector(test_folder + "b.txt");

    // correct equality constraints for normalization
    for (int i = 0; i < b_expected.size(); i++)
    {
        if (std::abs(b_expected(i) - Z->get_b()(i)) > 1e-3)
        {
            Ac_expected.row(i) *= Z->get_b()(i) / b_expected(i);
            Ab_expected.row(i) *= Z->get_b()(i) / b_expected(i);
            b_expected(i) = Z->get_b()(i);
        }
    }

    const HybZono Z_expected(Gc_expected, Gb_expected, c_expected, Ac_expected, Ab_expected, b_expected);

    // check that all matrices / vectors are approximately equal
    bool passed = true;
    passed &= Z->get_Gc().isApprox(Z_expected.get_Gc());
    passed &= Z->get_Gb().isApprox(Z_expected.get_Gb());
    passed &= Z->get_c().isApprox(Z_expected.get_c());
    passed &= Z->get_Ac().isApprox(Z_expected.get_Ac());
    passed &= Z->get_Ab().isApprox(Z_expected.get_Ab());
    passed &= Z->get_b().isApprox(Z_expected.get_b());

    test_assert(passed, "Minkowski sum test failed: computed result does not match expected result");

    return 0;
}