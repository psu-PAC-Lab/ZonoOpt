#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <filesystem>

using namespace ZonoOpt;

static std::filesystem::path tmp_dir()
{
    return std::filesystem::temp_directory_path() / "zonoopt_tests";
}

static void setup_tmp_dir()
{
    std::filesystem::create_directories(tmp_dir());
}

static void cleanup_tmp_dir()
{
    std::filesystem::remove_all(tmp_dir());
}

void test_zono()
{
    const ZonoPtr Z = make_regular_zono_2D(3., 12);
    const std::string filename = (tmp_dir() / "test_zono.json").string();
    to_json(*Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "test_zono: sets are not equal, expected " << *Z << ", got " << *Z_read;
    test_assert(hz_eq(*Z, *Z_read), ss.str());
    test_assert(Z_read->is_zono(), "test_zono: expected set to be of type Zono");
}

void test_hybzono()
{
    // hybrid zonotope
    Eigen::Matrix<zono_float, 2, 20> Gc;
    Gc << 0., 0.0859281, 0., 0., 0., 0., 0., 0.284171, 0.308149, 0., 0.116677, 0., 0., 0., 0., 0., 0., 0., 0., 0.798126,
         0.211588, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.511238, 0., 0.644165, 0., 0., 0., 0.28926, 0., 0., 0.;
    Eigen::Matrix<zono_float, 2, 5> Gb;
    Gb << 0, 0.235959, 0.26797, 0.669308, 0.757279,
            0,        0,        0,        0,        0;
    Eigen::Vector<zono_float, 2> c;
    c << 0.209747, 0.0100703;
    Eigen::Matrix<zono_float, 5, 20> Ac;
    Ac << 0,  0.731125,         0,  0.853555,   0.63719,  0.174854,         0,         0,  0.582327,         0,         0,         0,         0,         0,  0.575434, 0.0713598,         0,         0,         0,         0,
        0,         0,         0,   0.99319,   0.72748,         0,         0,         0,         0,         0,  0.355254,         0,  0.954118,         0,         0,         0,         0,         0,         0,         0,
        0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,  0.771462,         0,         0,         0,         0,         0,         0,         0,
  0.81856,         0,         0,  0.576755,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
        0,         0,         0,         0, 0.0988853,         0,         0,         0,  0.306937,  0.262899,         0,         0,         0,         0,         0,         0,         0,         0,   0.76945,         0;
    Eigen::Matrix<zono_float, 5, 5> Ab;
    Ab << 0,        0, 0.660926,        0,        0,
        0,        0,        0,        0,        0,
        0.742336, 0,        0, 0.474032,        0,
        0,        0,        0,        0, 0.560643,
        0.327432, 0,        0,        0,        0;
    Eigen::Vector<zono_float, 5> b;
    b << 0.243035,
         0.292617,
         0.610422,
         0.173898,
         0.702892;

    HybZono Z (Gc.sparseView(), Gb.sparseView(), c, Ac.sparseView(), Ab.sparseView(), b);

    // test json
    const std::string filename = (tmp_dir() / "test_hybzono.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "test_hybzono: sets are not equal, expected " << Z << ", got " << *Z_read;
    test_assert(hz_eq(Z, *Z_read), ss.str());
    test_assert(Z_read->is_hybzono(), "test_hybzono: expected set to be of type HybZono");
}

void test_conzono()
{
    // constrained zonotope
    Eigen::Matrix<zono_float, 1, 4> G;
    G << 3., 1., 0., 0.;
    Eigen::Vector<zono_float, 1> c;
    c << 8.;
    Eigen::Matrix<zono_float, 2, 4> A;
    A << 0.5, 0.1, 1., 0.,
         0., 0.5, 0., 0.5;
    Eigen::Vector<zono_float, 2> b;
    b << -0.5, -1.;

    ConZono Z (G.sparseView(), c, A.sparseView(), b);

    // test json
    const std::string filename = (tmp_dir() / "test_conzono.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "test_conzono: sets are not equal, expected " << Z << ", got " << *Z_read;
    test_assert(hz_eq(Z, *Z_read), ss.str());
    test_assert(Z_read->is_conzono(), "test_conzono: expected set to be of type ConZono");
}

void test_point()
{
    Eigen::Vector<zono_float, 4> c;
    c << 1., 2., 3., 4.;
    Point Z (c);

    // test json
    const std::string filename = (tmp_dir() / "test_point.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "test_point: sets are not equal, expected " << Z << ", got " << *Z_read;
    test_assert(hz_eq(Z, *Z_read), ss.str());
    test_assert(Z_read->is_point(), "test_point: expected set to be of type Point");
}

void test_empty_set()
{
    EmptySet Z (6);

    // test json
    const std::string filename = (tmp_dir() / "test_empty_set.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "test_empty_set: sets are not equal, expected " << Z << ", got " << *Z_read;
    test_assert(hz_eq(Z, *Z_read), ss.str());
    test_assert(Z_read->is_empty_set(), "test_empty_set: expected set to be of type EmptySet");
}

int main()
{
    setup_tmp_dir();

    test_hybzono();
    test_conzono();
    test_zono();
    test_point();
    test_empty_set();

    cleanup_tmp_dir();
    return 0;
}