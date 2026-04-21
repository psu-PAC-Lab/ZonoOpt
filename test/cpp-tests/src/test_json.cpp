#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <filesystem>

using namespace ZonoOpt;

static std::filesystem::path tmp_dir()
{
    return std::filesystem::temp_directory_path() / "zonoopt_tests";
}

class JsonTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        std::filesystem::create_directories(tmp_dir());
    }

    void TearDown() override
    {
        std::filesystem::remove_all(tmp_dir());
    }
};

TEST_F(JsonTest, Zono)
{
    const ZonoPtr Z = make_regular_zono_2D(3., 12);
    const std::string filename = (tmp_dir() / "test_zono.json").string();
    to_json(*Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "sets are not equal, expected " << *Z << ", got " << *Z_read;
    EXPECT_TRUE(hz_eq(*Z, *Z_read)) << ss.str();
    EXPECT_TRUE(Z_read->is_zono()) << "expected set to be of type Zono";
}

TEST_F(JsonTest, HybZono)
{
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

    const std::string filename = (tmp_dir() / "test_hybzono.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "sets are not equal, expected " << Z << ", got " << *Z_read;
    EXPECT_TRUE(hz_eq(Z, *Z_read)) << ss.str();
    EXPECT_TRUE(Z_read->is_hybzono()) << "expected set to be of type HybZono";
}

TEST_F(JsonTest, ConZono)
{
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

    const std::string filename = (tmp_dir() / "test_conzono.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "sets are not equal, expected " << Z << ", got " << *Z_read;
    EXPECT_TRUE(hz_eq(Z, *Z_read)) << ss.str();
    EXPECT_TRUE(Z_read->is_conzono()) << "expected set to be of type ConZono";
}

TEST_F(JsonTest, Point)
{
    Eigen::Vector<zono_float, 4> c;
    c << 1., 2., 3., 4.;
    Point Z (c);

    const std::string filename = (tmp_dir() / "test_point.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "sets are not equal, expected " << Z << ", got " << *Z_read;
    EXPECT_TRUE(hz_eq(Z, *Z_read)) << ss.str();
    EXPECT_TRUE(Z_read->is_point()) << "expected set to be of type Point";
}

TEST_F(JsonTest, EmptySet)
{
    EmptySet Z (6);

    const std::string filename = (tmp_dir() / "test_empty_set.json").string();
    to_json(Z, filename);

    const ZonoPtr Z_read = from_json(filename);

    std::stringstream ss;
    ss << "sets are not equal, expected " << Z << ", got " << *Z_read;
    EXPECT_TRUE(hz_eq(Z, *Z_read)) << ss.str();
    EXPECT_TRUE(Z_read->is_empty_set()) << "expected set to be of type EmptySet";
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
