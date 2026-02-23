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
    test_folder += "/support/";

    // load in conzono
    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "c.txt");
    const Eigen::SparseMatrix<zono_float>A = load_sparse_matrix(test_folder + "A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "b.txt");

    ConZono Z (G, c, A, b);

    // load direction and expected support value
    const Eigen::Vector<zono_float, -1> d = load_vector(test_folder + "d.txt");
    const zono_float s_expected = load_vector(test_folder + "sup.txt")(0);

    // compute support
    const zono_float s = Z.support(d);

    // compare results
    test_assert(std::abs(s - s_expected)/std::abs(s_expected) < 5e-2, "Support value does not match expected value");

    return 0;
}