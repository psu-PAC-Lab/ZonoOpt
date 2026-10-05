#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cstdlib>
#include <cmath>

using namespace ZonoOpt;

#define n_dims 3

zono_float f(const Eigen::Vector<double, n_dims>& x)
{
    return 2*std::pow(std::tan(x[0]), -2) + std::cos(x[1]/x[0])/3. + std::sin(x[0] + std::atan(x[2]))*std::sinh(x[0]) + std::exp(std::acosh(std::abs(x[1]) + 1)) - std::acos(x[0])*std::asin(x[1])/std::log(std::pow(x[2], 2));
}

Interval f_int(const Box& x)
{
    if (x.size() != n_dims)
        throw std::invalid_argument("Box must have size " + std::to_string(n_dims));

    Interval x0 = x.get_element(0);
    Interval x1 = x.get_element(1);
    Interval x2 = x.get_element(2);

    return 2.*(x0.tan().pow(-2)) + (x1/x0).cos()/3. + (x0 + x2.arctan()).sin()*x0.sinh() + (1. + x1.abs()).arccosh().exp() - (x0.arccos()*x1.arcsin())/(x2.pow(2)).log();
}

static void run_interval_test(zono_float x_min, zono_float x_max)
{
    Eigen::Vector<zono_float, -1> x_lb = Eigen::Vector<zono_float, n_dims>::Constant(x_min);
    Eigen::Vector<zono_float, -1> x_ub = Eigen::Vector<zono_float, n_dims>::Constant(x_max);
    const Box x (x_lb, x_ub);

    const Interval f_interval = f_int(x);

    srand(0);
    const int n_samples = 10000;

    for (int i=0; i < n_samples; ++i)
    {
        Eigen::Vector<double, n_dims> x_sample = Eigen::Vector<double, n_dims>::Random()*(x_max - x_min)/2. + Eigen::Vector<double, n_dims>::Constant((x_max + x_min)/2.);
        const zono_float f_sample = f(x_sample);

        std::stringstream err_ss;
        err_ss << "f(" << x_sample.transpose() << ") = " << f_sample << " not in " << f_interval;
        EXPECT_TRUE(f_interval.contains(f_sample)) << err_ss.str();
    }
}

TEST(IntervalArithmetic, PositiveRange)
{
    run_interval_test(0.1, 0.2);
}

TEST(IntervalArithmetic, NegativeRange)
{
    run_interval_test(-0.2, -0.001);
}

TEST(IntervalArithmetic, SpanningZero)
{
    run_interval_test(-1., 1.);
}

TEST(IntervalArithmetic, FractionalExponent)
{
    Interval a (0.5, 3.);
    Interval b = a.pow(456./123.);

    EXPECT_FALSE(b.is_empty()) << "test_exponent did not succeed";
    EXPECT_NEAR(b.lower(), std::pow(0.5, 456./123.), 1e-6) << "test_exponent lower bound incorrect";
    EXPECT_NEAR(b.upper(), std::pow(3., 456./123.), 1e-6) << "test_exponent upper bound incorrect";

    a = Interval(-3., -0.5);
    EXPECT_THROW(a.pow(456./123.), std::domain_error);
}

TEST(IntervalMatrix, FromTriplets)
{
    // entries at the same position are summed; positions without triplets are [0, 0]
    const std::vector<Eigen::Triplet<Interval>> triplets = {
        {0, 1, Interval(1, 2)}, {1, 2, Interval(-3, -1)}, {0, 1, Interval(0.5, 0.5)}};
    const IntervalMatrix M(2, 3, triplets);
    const auto vals = M.to_array();

    ASSERT_EQ(vals.size(), 2u);
    ASSERT_EQ(vals[0].size(), 3u);
    EXPECT_DOUBLE_EQ(vals[0][1].lower(), 1.5);
    EXPECT_DOUBLE_EQ(vals[0][1].upper(), 2.5);
    EXPECT_DOUBLE_EQ(vals[1][2].lower(), -3);
    EXPECT_DOUBLE_EQ(vals[1][2].upper(), -1);
    EXPECT_DOUBLE_EQ(vals[1][0].lower(), 0);
    EXPECT_DOUBLE_EQ(vals[1][0].upper(), 0);
}

TEST(IntervalMatrix, FromTripletsOutOfRangeThrows)
{
    const Interval iv(0, 1);
    for (const auto& [row, col] : std::vector<std::pair<int, int>>{{2, 0}, {0, 3}, {-1, 0}, {0, -1}})
    {
        const std::vector<Eigen::Triplet<Interval>> triplets = {{row, col, iv}};
        EXPECT_THROW(IntervalMatrix(2, 3, triplets), std::out_of_range) << "index (" << row << ", " << col << ")";
    }
}

TEST(IntervalArithmetic, FractionalExponentLargeDenominator)
{
    // exponents whose rational approximations need very large denominators
    const Interval a (0.5, 3.);
    constexpr double pi = 3.141592653589793;
    for (const zono_float p : {static_cast<zono_float>(pi), static_cast<zono_float>(0.50000000001)})
    {
        const Interval b = a.pow(p);
        EXPECT_NEAR(b.lower(), std::pow(0.5, p), 1e-6) << "lower bound incorrect for exponent " << p;
        EXPECT_NEAR(b.upper(), std::pow(3., p), 1e-6) << "upper bound incorrect for exponent " << p;
    }
}

TEST(IntervalArithmetic, Containment)
{
    const Interval a (0, 10);
    const Interval b (1, 2);
    EXPECT_TRUE(a.contains_set(b));
    EXPECT_FALSE(b.contains_set(a));
    EXPECT_TRUE(b <= a);
    EXPECT_TRUE(a >= b);
    EXPECT_TRUE(b == Interval(1, 2));
    EXPECT_FALSE(a == b);
}

TEST(IntervalArithmetic, EmptySetContainment)
{
    // the empty set is a subset of every set, and only the empty set is a subset of it
    const Interval e = Interval(0, 1).intersect(Interval(2, 3));
    const Interval a (0, 10);
    ASSERT_TRUE(e.is_empty());
    EXPECT_TRUE(e == e);
    EXPECT_TRUE(a.contains_set(e));
    EXPECT_TRUE(e <= a);
    EXPECT_FALSE(e.contains_set(a));
    EXPECT_FALSE(e == a);

    // width of an interval matrix with an empty element is not a number
    const IntervalMatrix M(1, 2, {{0, 0, e}, {0, 1, a}});
    EXPECT_TRUE(std::isnan(M.width()));
}

class IntervalMatrixScalarOps : public ::testing::Test
{
protected:
    // 2 x 2 matrix with only element (0, 0) stored, the others are implicit zeros
    IntervalMatrix M = IntervalMatrix(2, 2, {{0, 0, Interval(1, 2)}});

    static void expect_element(const IntervalMatrix& A, const int i, const int j, const Interval& expected,
                               const std::string& name)
    {
        const Interval val = A.to_array()[i][j];
        EXPECT_NEAR(val.lower(), expected.lower(), 1e-12) << name << " (" << i << ", " << j << ")";
        EXPECT_NEAR(val.upper(), expected.upper(), 1e-12) << name << " (" << i << ", " << j << ")";
    }
};

TEST_F(IntervalMatrixScalarOps, ZeroScalarLeavesImplicitZeros)
{
    for (const auto& [A, name] : std::vector<std::pair<IntervalMatrix, std::string>>{
             {M + 0., "M + 0"}, {M - Interval(0, 0), "M - [0, 0]"}, {M * 2., "M * 2"}})
    {
        expect_element(A, 0, 0, name == "M * 2" ? Interval(2, 4) : Interval(1, 2), name);
        expect_element(A, 0, 1, Interval(0, 0), name);
        expect_element(A, 1, 1, Interval(0, 0), name);
    }
}

TEST_F(IntervalMatrixScalarOps, ScalarOpsApplyToImplicitZeros)
{
    expect_element(M + 1., 1, 0, Interval(1, 1), "M + 1");
    expect_element(M + Interval(1, 2), 0, 1, Interval(1, 2), "M + [1, 2]");
    expect_element(M - 1., 1, 1, Interval(-1, -1), "M - 1");
    expect_element(M - Interval(1, 2), 1, 1, Interval(-2, -1), "M - [1, 2]");
    expect_element(3. - M, 1, 0, Interval(3, 3), "3 - M");
    expect_element(Interval(1, 2) - M, 0, 1, Interval(1, 2), "[1, 2] - M");
    expect_element(M + 1., 0, 0, Interval(2, 3), "M + 1");

    // division by an implicit zero matches division by [0, 0], not [0, 0]
    const Interval div_zero = 1. / Interval(0, 0);
    const Interval val = (1. / M).to_array()[0][1];
    EXPECT_EQ(val.is_empty(), div_zero.is_empty()) << "1 / M (0, 1)";
    if (!div_zero.is_empty())
    {
        EXPECT_EQ(val.lower(), div_zero.lower()) << "1 / M (0, 1)";
        EXPECT_EQ(val.upper(), div_zero.upper()) << "1 / M (0, 1)";
    }
    EXPECT_FALSE(val.lower() == 0 && val.upper() == 0) << "1 / M (0, 1) must not be [0, 0]";
}

TEST(IntervalMatrix, FromDenseMatrices)
{
    // zero elements are implicit, others are stored
    Eigen::Matrix<Interval, -1, -1> mat(2, 2);
    mat << Interval(1, 2), Interval(0, 0), Interval(0, 0), Interval(-3, -1);
    Eigen::Matrix<zono_float, -1, -1> lb(2, 2), ub(2, 2);
    lb << 1, 0, 0, -3;
    ub << 2, 0, 0, -1;

    for (const auto& [M, name] : std::vector<std::pair<IntervalMatrix, std::string>>{
             {IntervalMatrix(mat), "Matrix<Interval>"}, {IntervalMatrix(lb, ub), "lower / upper bounds"}})
    {
        const auto vals = M.to_array();
        EXPECT_TRUE(vals[0][0] == Interval(1, 2)) << name;
        EXPECT_TRUE(vals[0][1] == Interval(0, 0)) << name;
        EXPECT_TRUE(vals[1][1] == Interval(-3, -1)) << name;
        EXPECT_FALSE(M.is_empty()) << name;
    }
}

TEST(IntervalMatrix, FromDenseMatricesKeepsEmptyElements)
{
    // empty elements have NaN bounds and must not be treated as zero
    const zono_float nan = std::numeric_limits<zono_float>::quiet_NaN();
    Eigen::Matrix<Interval, -1, -1> mat(1, 2);
    mat << Interval(0, 1).intersect(Interval(2, 3)), Interval(1, 2);
    Eigen::Matrix<zono_float, -1, -1> lb(1, 2), ub(1, 2);
    lb << nan, 1;
    ub << nan, 2;

    for (const auto& [M, name] : std::vector<std::pair<IntervalMatrix, std::string>>{
             {IntervalMatrix(mat), "Matrix<Interval>"}, {IntervalMatrix(lb, ub), "lower / upper bounds"}})
    {
        EXPECT_TRUE(M.is_empty()) << name;
        EXPECT_TRUE(M.to_array()[0][0].is_empty()) << name;
    }
}
