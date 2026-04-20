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
    Eigen::Vector<zono_float, -1> d = load_vector(test_folder + "d.txt");
    const zono_float s_expected = load_vector(test_folder + "sup.txt")(0);

    // compute support
    const zono_float s = Z.support(d);

    // compare results
    test_assert(std::abs(s - s_expected)/std::abs(s_expected) < 5e-2, "Support value does not match expected value");

    // test that Zono and ConZono return matching support values and factors
    d.resize(5);
    d << 1., -1., 0.5, -0.5, 1.5;
    OptSettings settings;
    settings.eps_prim = 1e-3;
    settings.eps_dual = 1e-3;
    settings.rho = 1.;
    for (unsigned int i=0; i<10; ++i)
    {
        std::string filename = test_folder + "zono_" + std::to_string(i) + ".json";
        auto Zjson = from_json(filename);
        Zono* Zz_tmp = dynamic_cast<Zono*>(Zjson.get());
        if (Zz_tmp == nullptr)
        {
            std::cerr << "Error: loaded set from " << filename << " is not a Zono." << std::endl;
            return 1;
        }
        Zjson.release();
        std::unique_ptr<Zono> Zz (Zz_tmp); 

        ConZono Zc(Zz->get_G(), Zz->get_c(), Zz->get_A(), Zz->get_b());

        auto sol1 = std::make_shared<OptSolution>();
        auto sol2 = std::make_shared<OptSolution>();
        const zono_float s1 = Zz->support(d, settings, &sol1);
        const zono_float s2 = Zc.support(d, settings, &sol2);

        std::stringstream ss;

        ss << "Zono and ConZono support values do not match" << "\nZono support: " << s1 << "\nConZono support: " << s2;
        ss << "Z2: " << *Zz << ", Zc: " << Zc;
        ss << "sol1: " << sol1->print() << "\nsol2: " << sol2->print();
        test_assert(std::abs(s1 - s2) < 1e-2, ss.str());

        ss.str("");
        ss << "Zono and ConZono solutions do not produce same support value when applied to generators and center" << 
            "\nZono support: " << d.dot(Zz->get_G() * sol1->z + Zz->get_c()) << "\nConZono support: " << 
            d.dot(Zc.get_G() * sol2->z + Zc.get_c());
        test_assert(std::abs(d.dot(Zz->get_G() * sol1->z + Zz->get_c()) - d.dot(Zc.get_G() * sol2->z + Zc.get_c())) < 1e-2, ss.str());
    }

    return 0;
}