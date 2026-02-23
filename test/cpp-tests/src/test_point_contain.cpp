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
    test_folder += "/point_contain/";

    // load in conzono
    const Eigen::SparseMatrix<zono_float> G = load_sparse_matrix(test_folder + "G.txt");
    const Eigen::Vector<zono_float, -1> c = load_vector(test_folder + "c.txt");
    const Eigen::SparseMatrix<zono_float>A = load_sparse_matrix(test_folder + "A.txt");
    const Eigen::Vector<zono_float, -1> b = load_vector(test_folder + "b.txt");

    ConZono Z (G, c, A, b);

    // load point in set
    const Eigen::Vector<zono_float, -1> x_c = load_vector(test_folder + "x_c.txt");

    // load point not in set
    const Eigen::Vector<zono_float, -1> x_n = load_vector(test_folder + "x_n.txt");

    // check correct classification of containment
    test_assert(Z.contains_point(x_c), "Expected Z to contain x_c");
    test_assert(!Z.contains_point(x_n), "Expected Z to not contain x_n");

    return 0;
}