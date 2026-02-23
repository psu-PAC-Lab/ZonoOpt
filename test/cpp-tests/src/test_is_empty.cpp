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
    test_folder += "/is_empty/";

    // load in feasible conzono
    Eigen::SparseMatrix<zono_float> G, A;
    Eigen::Vector<zono_float, -1> c, b;
    G = load_sparse_matrix(test_folder + "f_G.txt");
    c = load_vector(test_folder + "f_c.txt");
    A = load_sparse_matrix(test_folder + "f_A.txt");
    b = load_vector(test_folder + "f_b.txt");

    const ConZono Zf (G, c, A, b);

    // load in infeasible conzono
    G = load_sparse_matrix(test_folder + "i_G.txt");
    c = load_vector(test_folder + "i_c.txt");
    A = load_sparse_matrix(test_folder + "i_A.txt");
    b = load_vector(test_folder + "i_b.txt");

    const ConZono Zi (G, c, A, b);

    // check if empty
    test_assert(!Zf.is_empty(), "Expected Zf to be non-empty");
    test_assert(Zi.is_empty(), "Expected Zi to be empty");

    return 0;
}