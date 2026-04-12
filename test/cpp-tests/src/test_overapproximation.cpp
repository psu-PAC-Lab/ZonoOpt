#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

void test_reduce_order(const std::string& filename)
{
    // from file
    ZonoPtr Z = from_json(filename);

    // convert to Zono
    if (!Z->is_zono())
    {
        throw std::invalid_argument("Expected set to be a zonotope for reduce order test.");
    }
    
    std::unique_ptr<Zono> Z_zono (dynamic_cast<Zono*>(Z.release()));

    // project point onto set
    Eigen::Vector<zono_float, -1> p (Z_zono->get_n());
    p.setZero();
    const auto p_proj = Z_zono->project_point(p);

    // reduce order
    const auto Zr = Z_zono->reduce_order(10);

    // make sure point is contained in reduced set
    test_assert(Zr->contains_point(p_proj), "Reduced zonotope does not contain projected point.");
}

void test_constraint_reduction(const std::string& filename)
{
    // from file
    ZonoPtr Z = from_json(filename);

    // convert to ConZono
    if (!Z->is_conzono())
    {
        throw std::invalid_argument("Expected set to be a constrained zonotope for constraint reduction test.");
    }
    
    std::unique_ptr<ConZono> Z_conzono (dynamic_cast<ConZono*>(Z.release()));

    // project point onto set
    Eigen::Vector<zono_float, -1> p (Z_conzono->get_n());
    p.setZero();
    const auto p_proj = Z_conzono->project_point(p);

    // constraint reduction
    const auto Z_cr = Z_conzono->constraint_reduction();

    // make sure point is contained in reduced set
    test_assert(Z_cr->contains_point(p_proj), "Reduced constrained zonotope does not contain projected point.");
}

void test_to_zono_approx(const std::string& filename)
{
    // from file
    ZonoPtr Z = from_json(filename);

    // convert to ConZono
    if (!Z->is_conzono())
    {
        throw std::invalid_argument("Expected set to be a constrained zonotope for zonotope approximation test.");
    }
    
    std::unique_ptr<ConZono> Z_conzono (dynamic_cast<ConZono*>(Z.release()));

    // project point onto set
    Eigen::Vector<zono_float, -1> p (Z_conzono->get_n());
    p.setZero();
    const auto p_proj = Z_conzono->project_point(p);

    // zonotope approximation
    const auto Z_approx = Z_conzono->to_zono_approx();

    // make sure point is contained in approximated set
    test_assert(Z_approx->contains_point(p_proj), "Zonotope approximation does not contain projected point.");
}

void test_convex_relaxation(const std::string& filename)
{
    // from file
    ZonoPtr Z = from_json(filename);

    // project point onto set
    Eigen::Vector<zono_float, -1> p (Z->get_n());
    p.setZero();
    const auto p_proj = Z->project_point(p);

    // convex relaxation
    const auto Z_relax = Z->convex_relaxation();

    // make sure point is contained in relaxed set
    test_assert(Z_relax->contains_point(p_proj), "Convex relaxation does not contain projected point.");
}


int main(int argc, char* argv[])
{
    // input: directory where unit test data resides
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <test data folder>" << std::endl;
        return 1;
    }
    std::string test_folder = argv[1];
    test_folder += "/overapproximation/";

    // run tests
    test_reduce_order(test_folder + "rand_zono.json");
    test_constraint_reduction(test_folder + "rand_conzono.json");
    test_to_zono_approx(test_folder + "rand_conzono.json");
    test_convex_relaxation(test_folder + "rand_hybzono.json");

    return 0;
}