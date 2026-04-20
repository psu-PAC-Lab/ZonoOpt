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

    // test that Zono and ConZono return matching support values and factors
    std::mt19937 rand_gen(42);
    OptSettings settings;
    settings.eps_prim = 1e-3;
    settings.eps_dual = 1e-3;
    settings.rho = 1.;

    Zono Z2 = random_zono(5, 10, 0.5, -1., 1., rand_gen);
    ConZono Zc(Z2.get_G(), Z2.get_c(), Z2.get_A(), Z2.get_b());
    const Eigen::Vector<zono_float, -1> d2 = random_vector(5, -1., 1., rand_gen);

    auto sol1 = std::make_shared<OptSolution>();
    auto sol2 = std::make_shared<OptSolution>();
    const zono_float s1 = Z2.support(d2, settings, &sol1);
    const zono_float s2 = Zc.support(d2, settings, &sol2);

    test_assert(std::abs(s1 - s2) < 1e-2, "Zono and ConZono support values do not match");
    test_assert((Z2.get_G() * sol1->z - Zc.get_G() * sol2->z).norm() < 1e-2, "Zono and ConZono solution factors do not match");

    return 0;
}