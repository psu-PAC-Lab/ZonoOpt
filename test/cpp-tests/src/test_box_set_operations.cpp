#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace ZonoOpt;

namespace
{
    // enumerates all 2^n corner points of a box (n assumed small, <=4, in these tests)
    std::vector<Eigen::Vector<zono_float, -1>> box_vertices(const Box& b)
    {
        const int n = static_cast<int>(b.size());
        std::vector<Eigen::Vector<zono_float, -1>> verts;
        const int n_verts = 1 << n;
        for (int mask = 0; mask < n_verts; ++mask)
        {
            Eigen::Vector<zono_float, -1> v(n);
            for (int i = 0; i < n; ++i)
                v(i) = (mask & (1 << i)) ? b.upper()(i) : b.lower()(i);
            verts.push_back(v);
        }
        return verts;
    }
}

TEST(BoxSetOperations, IsEmptyDetectsBoostEmpty)
{
    // intersection of disjoint boxes: boost encodes the resulting interval as [NaN, NaN],
    // which must be detected as empty rather than compared numerically
    Box b1(Eigen::Vector<zono_float, -1>::Constant(1, 0.), Eigen::Vector<zono_float, -1>::Constant(1, 1.));
    Box b2(Eigen::Vector<zono_float, -1>::Constant(1, 2.), Eigen::Vector<zono_float, -1>::Constant(1, 3.));
    EXPECT_TRUE(b1.intersect(b2).is_empty());

    // raw lb > ub encoding must also be detected
    Box b3(Eigen::Vector<zono_float, -1>::Constant(1, 1.), Eigen::Vector<zono_float, -1>::Constant(1, 0.));
    EXPECT_TRUE(b3.is_empty());

    // a normal box is not empty
    EXPECT_FALSE(b1.is_empty());

    // zero-dimension box is empty by convention
    EXPECT_TRUE(Box(0).is_empty());
}

TEST(BoxSetOperations, CartesianProduct)
{
    Eigen::Vector<zono_float, -1> lb1(1), ub1(1), lb2(2), ub2(2);
    lb1 << 0.; ub1 << 1.;
    lb2 << 2., 4.; ub2 << 3., 5.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    Box prod = b1.cartesian_product(b2);
    EXPECT_EQ(prod.size(), 3u);

    // definitional: (x, y) is in b1 x b2 for every x in b1, y in b2
    for (const auto& x : box_vertices(b1))
    {
        for (const auto& y : box_vertices(b2))
        {
            Eigen::Vector<zono_float, -1> xy(3);
            xy << x, y;
            EXPECT_TRUE(prod.contains(xy));
        }
    }

    // definitional converse: a point violating either factor is excluded
    Eigen::Vector<zono_float, -1> outside(3);
    outside << -1., 2., 4.; // first component not in b1 = [0,1]
    EXPECT_FALSE(prod.contains(outside));

    EXPECT_TRUE(cartesian_product(b1, b2) == prod);
}

TEST(BoxSetOperations, CartesianProductEmptyPropagates)
{
    Box b1(Eigen::Vector<zono_float, -1>::Constant(1, 0.), Eigen::Vector<zono_float, -1>::Constant(1, 1.));
    Box b_empty(Eigen::Vector<zono_float, -1>::Constant(1, 1.), Eigen::Vector<zono_float, -1>::Constant(1, 0.));

    EXPECT_TRUE(b1.cartesian_product(b_empty).is_empty());
    EXPECT_TRUE(b_empty.cartesian_product(b1).is_empty());
}

TEST(BoxSetOperations, MinkowskiSum)
{
    Eigen::Vector<zono_float, -1> lb1(2), ub1(2), lb2(2), ub2(2);
    lb1 << 0., -1.; ub1 << 10., 4.;
    lb2 << 1., -2.; ub2 << 2., 3.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    Box sum = b1.minkowski_sum(b2);
    EXPECT_TRUE(sum == (b1 + b2));
    EXPECT_TRUE(sum == minkowski_sum(b1, b2));

    // definitional: x + y is in the sum for every x in b1, y in b2
    for (const auto& x : box_vertices(b1))
        for (const auto& y : box_vertices(b2))
            EXPECT_TRUE(sum.contains(x + y));

    // cross-check against the zonotope path
    auto Z1 = interval_2_zono(b1);
    auto Z2 = interval_2_zono(b2);
    const auto Z_sum = minkowski_sum(*Z1, *Z2);
    const Box bb = Z_sum->bounding_box();
    EXPECT_TRUE(vector_eq(bb.lower(), sum.lower()));
    EXPECT_TRUE(vector_eq(bb.upper(), sum.upper()));
}

TEST(BoxSetOperations, MinkowskiSumEmptyPropagates)
{
    // deliberate divergence from HybZono::minkowski_sum, which treats an empty
    // operand as the identity; Box::minkowski_sum stays a pure alias for operator+
    Box b1(Eigen::Vector<zono_float, -1>::Constant(1, 0.), Eigen::Vector<zono_float, -1>::Constant(1, 1.));
    Box b_empty(Eigen::Vector<zono_float, -1>::Constant(1, 1.), Eigen::Vector<zono_float, -1>::Constant(1, 0.));

    EXPECT_TRUE(b1.minkowski_sum(b_empty).is_empty());
    EXPECT_TRUE((b1 + b_empty).is_empty());
}

TEST(BoxSetOperations, PontryDiff)
{
    Eigen::Vector<zono_float, -1> lb1(2), ub1(2), lb2(2), ub2(2);
    lb1 << 0., -4.; ub1 << 10., 4.;
    lb2 << 1., -1.; ub2 << 2., 2.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    Box diff = b1.pontry_diff(b2);
    EXPECT_TRUE(diff.get_element(0) == Interval(-1., 8.));
    EXPECT_TRUE(diff.get_element(1) == Interval(-3., 2.));

    // definitional: (b1 (-) b2) + b2 is contained in b1, i.e., x + y in b1 for all
    // x in the difference and y in b2
    EXPECT_TRUE(b1.contains_set(diff.minkowski_sum(b2)));

    // must differ from operator- (interval subtraction gives [a-d, b-c], not the
    // Pontryagin difference [a-c, b-d])
    EXPECT_FALSE(diff == (b1 - b2));

    EXPECT_TRUE(pontry_diff(b1, b2) == diff);
}

TEST(BoxSetOperations, PontryDiffEmptyResult)
{
    // subtrahend wider than minuend => empty result
    Eigen::Vector<zono_float, -1> lb1(1), ub1(1), lb2(1), ub2(1);
    lb1 << 0.; ub1 << 1.;
    lb2 << -5.; ub2 << 5.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    EXPECT_TRUE(b1.pontry_diff(b2).is_empty());
}

TEST(BoxSetOperations, PontryDiffEmptyInputReturnsMinuend)
{
    // matches HybZono::pontry_diff: either input empty => result equals the minuend
    Box b1(Eigen::Vector<zono_float, -1>::Constant(1, 0.), Eigen::Vector<zono_float, -1>::Constant(1, 1.));
    Box b_empty(Eigen::Vector<zono_float, -1>::Constant(1, 1.), Eigen::Vector<zono_float, -1>::Constant(1, 0.));

    EXPECT_TRUE(b1.pontry_diff(b_empty) == b1);
    EXPECT_TRUE(b_empty.pontry_diff(b1).is_empty());
}

TEST(BoxSetOperations, PontryDiffCrossCheckZono)
{
    Eigen::Vector<zono_float, -1> lb1(2), ub1(2), lb2(2), ub2(2);
    lb1 << 0., -1.; ub1 << 10., 4.;
    lb2 << 1., -2.; ub2 << 2., 3.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    auto Z1 = interval_2_zono(b1);
    auto Z2 = interval_2_zono(b2);
    const auto Z_diff = pontry_diff(*Z1, *Z2, true);
    const Box bb = Z_diff->bounding_box();

    const Box diff = b1.pontry_diff(b2);
    for (int i = 0; i < static_cast<int>(diff.size()); ++i)
    {
        EXPECT_NEAR(bb.lower()(i), diff.lower()(i), 1e-2);
        EXPECT_NEAR(bb.upper()(i), diff.upper()(i), 1e-2);
    }
}

TEST(BoxSetOperations, ProjectOntoDims)
{
    Eigen::Vector<zono_float, -1> lb(3), ub(3);
    lb << 0., 1., 2.;
    ub << 10., 11., 12.;
    Box b(lb, ub);

    const std::vector<int> dims{0, 2};
    Box proj = b.project_onto_dims(dims);
    EXPECT_EQ(proj.size(), 2u);

    // definitional: selecting components dims of any point in b lands in the
    // projection, and every point of the projection is realized by some point in b
    for (const auto& v : box_vertices(b))
    {
        Eigen::Vector<zono_float, -1> selected(2);
        selected << v(0), v(2);
        EXPECT_TRUE(proj.contains(selected));
    }
    for (const auto& w : box_vertices(proj))
    {
        Eigen::Vector<zono_float, -1> extended(3);
        extended << w(0), b.center()(1), w(1);
        EXPECT_TRUE(b.contains(extended));
    }

    // reordered
    Box reordered = b.project_onto_dims({2, 0});
    EXPECT_TRUE(reordered.get_element(0) == Interval(2., 12.));
    EXPECT_TRUE(reordered.get_element(1) == Interval(0., 10.));

    // repeated
    Box repeated = b.project_onto_dims({1, 1});
    EXPECT_TRUE(repeated.get_element(0) == Interval(1., 11.));
    EXPECT_TRUE(repeated.get_element(1) == Interval(1., 11.));

    EXPECT_TRUE(project_onto_dims(b, dims) == proj);
}

TEST(BoxSetOperations, ProjectOntoDimsThrows)
{
    Box b(Eigen::Vector<zono_float, -1>::Zero(3), Eigen::Vector<zono_float, -1>::Ones(3));
    EXPECT_THROW(b.project_onto_dims({-1}), std::invalid_argument);
    EXPECT_THROW(b.project_onto_dims({3}), std::invalid_argument);
}

TEST(BoxSetOperations, ProjectOntoDimsEmptyInput)
{
    Eigen::Vector<zono_float, -1> lb(2), ub(2);
    lb << 0., 5.;
    ub << 1., 4.; // second component empty
    Box b(lb, ub);
    EXPECT_TRUE(b.is_empty());

    // projecting away the empty component must not un-empty the box
    EXPECT_TRUE(b.project_onto_dims({0}).is_empty());
}

TEST(BoxSetOperations, AffineMapOverloadsAgree)
{
    Eigen::Vector<zono_float, -1> lb(2), ub(2);
    lb << 0., -1.;
    ub << 10., 4.;
    Box b(lb, ub);

    Eigen::Matrix<zono_float, -1, -1> A_dense(2, 2);
    A_dense << 1., 2., -1., 0.5;
    Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_csr = A_dense.sparseView();
    Eigen::SparseMatrix<zono_float> A_csc = A_dense.sparseView();

    Eigen::Vector<zono_float, -1> s(2);
    s << 1., -2.;

    Box out_dense = b.affine_map(A_dense, s);
    Box out_csr = b.affine_map(A_csr, s);
    Box out_csc = b.affine_map(A_csc, s);

    EXPECT_TRUE(vector_eq(out_dense.lower(), out_csr.lower()));
    EXPECT_TRUE(vector_eq(out_dense.upper(), out_csr.upper()));
    EXPECT_TRUE(vector_eq(out_dense.lower(), out_csc.lower()));
    EXPECT_TRUE(vector_eq(out_dense.upper(), out_csc.upper()));

    // definitional: R*x + s is in the affine map's output for every x in the box
    for (const auto& x : box_vertices(b))
        EXPECT_TRUE(out_dense.contains(A_dense * x + s));

    // affine_map(A, s) == linear_map(A) + s
    Box expected = b.linear_map(A_dense) + s;
    EXPECT_TRUE(vector_eq(out_dense.lower(), expected.lower()));
    EXPECT_TRUE(vector_eq(out_dense.upper(), expected.upper()));

    // empty s defaults to zero offset
    Box no_offset = b.affine_map(A_dense);
    EXPECT_TRUE(vector_eq(no_offset.lower(), b.linear_map(A_dense).lower()));

    // free function (ColMajor) agrees with member
    Box free_result = affine_map(b, A_csc, s);
    EXPECT_TRUE(vector_eq(free_result.lower(), out_csc.lower()));
}

TEST(BoxSetOperations, AffineMapEmptyInput)
{
    Eigen::Vector<zono_float, -1> lb(2), ub(2);
    lb << 0., 5.;
    ub << 1., 4.; // second component empty
    Box b(lb, ub);

    Eigen::Matrix<zono_float, -1, -1> A(1, 2);
    A << 1., 0.; // zero coefficient on the empty column
    EXPECT_TRUE(b.affine_map(A).is_empty());
}

TEST(BoxSetOperations, Support)
{
    Eigen::Vector<zono_float, -1> lb(3), ub(3);
    lb << -1., 2., -3.;
    ub << 4., 5., 1.;
    Box b(lb, ub);

    Eigen::Vector<zono_float, -1> d(3);
    d << 1., -1., 2.;

    // definitional: support(d) = max_{x in b} <x, d>; for a box this max is
    // attained at a vertex, so brute-force enumeration of vertices computes it exactly
    zono_float expected = -std::numeric_limits<zono_float>::infinity();
    for (const auto& v : box_vertices(b))
        expected = std::max(expected, static_cast<zono_float>(v.dot(d)));
    EXPECT_NEAR(b.support(d), expected, 1e-9);

    // cross-check against zonotope support
    auto Z = interval_2_zono(b);
    EXPECT_NEAR(b.support(d), Z->support(d), 1e-6);
}

TEST(BoxSetOperations, IntervalHullMany)
{
    Eigen::Vector<zono_float, -1> lb1(1), ub1(1), lb2(1), ub2(1), lb3(1), ub3(1);
    lb1 << 0.; ub1 << 1.;
    lb2 << -2.; ub2 << 0.5;
    lb3 << 3.; ub3 << 4.;
    Box b1(lb1, ub1), b2(lb2, ub2), b3(lb3, ub3);

    Box hull = interval_hull(std::vector<Box>{b1, b2, b3});

    // definitional: the hull contains every input box and is the tightest such box
    EXPECT_TRUE(hull.contains_set(b1));
    EXPECT_TRUE(hull.contains_set(b2));
    EXPECT_TRUE(hull.contains_set(b3));
    EXPECT_TRUE(hull.get_element(0) == Interval(-2., 4.));

    // single box
    Box single = interval_hull(std::vector<Box>{b1});
    EXPECT_TRUE(single.get_element(0) == b1.get_element(0));

    // empty vector throws
    EXPECT_THROW(interval_hull(std::vector<Box>{}), std::invalid_argument);

    // list containing empties: empties are ignored
    Box b_empty(Eigen::Vector<zono_float, -1>::Constant(1, 1.), Eigen::Vector<zono_float, -1>::Constant(1, 0.));
    Box hull_with_empty = interval_hull(std::vector<Box>{b_empty, b1, b2});
    EXPECT_TRUE(hull_with_empty.get_element(0) == interval_hull(std::vector<Box>{b1, b2}).get_element(0));

    // all-empty => empty
    EXPECT_TRUE(interval_hull(std::vector<Box>{b_empty, b_empty}).is_empty());

    // two-arg free function agrees with member
    EXPECT_TRUE(interval_hull(b1, b2) == b1.interval_hull(b2));
}

TEST(BoxSetOperations, Intersection)
{
    Eigen::Vector<zono_float, -1> lb1(2), ub1(2), lb2(2), ub2(2);
    lb1 << 0., 0.; ub1 << 10., 10.;
    lb2 << 5., -5.; ub2 << 15., 5.;
    Box b1(lb1, ub1), b2(lb2, ub2);

    Box result = intersection(b1, b2);
    EXPECT_TRUE(result == b1.intersect(b2));
    EXPECT_TRUE(result == (b1 & b2));
    EXPECT_TRUE(result.get_element(0) == Interval(5., 10.));
    EXPECT_TRUE(result.get_element(1) == Interval(0., 5.));

    // definitional: every vertex of the result is in both b1 and b2
    for (const auto& v : box_vertices(result))
    {
        EXPECT_TRUE(b1.contains(v));
        EXPECT_TRUE(b2.contains(v));
    }

    // definitional converse: a point in both b1 and b2 lies in the result
    Eigen::Vector<zono_float, -1> p(2);
    p << 7., 2.;
    ASSERT_TRUE(b1.contains(p) && b2.contains(p));
    EXPECT_TRUE(result.contains(p));

    // a point in b1 but not b2 is excluded from the result
    Eigen::Vector<zono_float, -1> q(2);
    q << 1., 1.;
    ASSERT_TRUE(b1.contains(q) && !b2.contains(q));
    EXPECT_FALSE(result.contains(q));
}

TEST(BoxSetOperations, IntersectionDisjoint)
{
    Box b1(Eigen::Vector<zono_float, -1>::Constant(1, 0.), Eigen::Vector<zono_float, -1>::Constant(1, 1.));
    Box b2(Eigen::Vector<zono_float, -1>::Constant(1, 2.), Eigen::Vector<zono_float, -1>::Constant(1, 3.));
    EXPECT_TRUE(intersection(b1, b2).is_empty());
}

TEST(BoxSetOperations, IntersectionEmptyInput)
{
    Box b1(Eigen::Vector<zono_float, -1>::Constant(2, 0.), Eigen::Vector<zono_float, -1>::Constant(2, 1.));
    Box b_empty(Eigen::Vector<zono_float, -1>::Constant(2, 1.), Eigen::Vector<zono_float, -1>::Constant(2, 0.));

    Box result = intersection(b1, b_empty);
    EXPECT_TRUE(result.is_empty());
    EXPECT_EQ(result.size(), b1.size());

    // also with an explicit (identity) R
    Eigen::SparseMatrix<zono_float> R(2, 2);
    R.setIdentity();
    Box result_R = intersection(b1, b_empty, R);
    EXPECT_TRUE(result_R.is_empty());
    EXPECT_EQ(result_R.size(), b1.size());
}

TEST(BoxSetOperations, IntersectionThrows)
{
    Box b1(Eigen::Vector<zono_float, -1>::Zero(2), Eigen::Vector<zono_float, -1>::Ones(2));
    Box b2(Eigen::Vector<zono_float, -1>::Zero(3), Eigen::Vector<zono_float, -1>::Ones(3));

    // default (identity) R requires matching dimensions
    EXPECT_THROW(intersection(b1, b2), std::invalid_argument);

    // explicit R with the wrong shape (should be Z2.size() x Z1.size() = 3x2)
    Eigen::SparseMatrix<zono_float> R(3, 3);
    R.setIdentity();
    EXPECT_THROW(intersection(b1, b2, R), std::invalid_argument);
}

TEST(BoxSetOperations, IntersectionSelectionRMatchesOverDims)
{
    Eigen::Vector<zono_float, -1> lb1(3), ub1(3);
    lb1 << 0., 1., 2.; ub1 << 10., 11., 12.;
    Box b1(lb1, ub1);

    Eigen::Vector<zono_float, -1> lb2(2), ub2(2);
    lb2 << 5., 3.; ub2 << 9., 20.;
    Box b2(lb2, ub2);

    const std::vector<int> dims{0, 2};

    // selection matrix R (2x3): row k picks out dims[k]
    Eigen::SparseMatrix<zono_float> R(static_cast<Eigen::Index>(dims.size()), static_cast<Eigen::Index>(b1.size()));
    std::vector<Eigen::Triplet<zono_float>> trips;
    for (size_t k = 0; k < dims.size(); ++k)
        trips.emplace_back(static_cast<int>(k), dims[k], static_cast<zono_float>(1.));
    R.setFromTriplets(trips.begin(), trips.end());

    // a selection R makes the generalized intersection exact, so it must agree exactly
    // with intersection_over_dims (both member and free forms)
    Box via_R = intersection(b1, b2, R);
    Box via_member_over_dims = b1.intersection_over_dims(b2, dims);
    Box via_free_over_dims = intersection_over_dims(b1, b2, dims);

    EXPECT_TRUE(via_R == via_member_over_dims);
    EXPECT_TRUE(via_member_over_dims == via_free_over_dims);
}

TEST(BoxSetOperations, IntersectionGeneralR)
{
    // R is not a selection matrix: constrain Z1 (2D) via x0 + x1 in Z2 (1D)
    Eigen::Vector<zono_float, -1> lb1(2), ub1(2);
    lb1 << 0., 0.; ub1 << 10., 10.;
    Box b1(lb1, ub1);

    Box b2(Eigen::Vector<zono_float, -1>::Constant(1, 5.), Eigen::Vector<zono_float, -1>::Constant(1, 6.));

    Eigen::SparseMatrix<zono_float> R(1, 2);
    std::vector<Eigen::Triplet<zono_float>> trips{{0, 0, static_cast<zono_float>(1.)}, {0, 1, static_cast<zono_float>(1.)}};
    R.setFromTriplets(trips.begin(), trips.end());

    Box result = intersection(b1, b2, R);

    // soundness: the contractor only narrows, so the result must stay within b1
    EXPECT_TRUE(b1.contains_set(result));

    // soundness: every point of b1 satisfying x0 + x1 in b2 must remain in the result,
    // since the interval contractor is guaranteed not to remove valid points
    for (zono_float x0 = 0.; x0 <= 10.; x0 += 1.)
    {
        for (zono_float x1 = 0.; x1 <= 10.; x1 += 1.)
        {
            if (x0 + x1 < 5. || x0 + x1 > 6.)
                continue;
            Eigen::Vector<zono_float, -1> x(2);
            x << x0, x1;
            EXPECT_TRUE(result.contains(x));
        }
    }

    // tightness cross-check: the exact generalized intersection on the zonotope path is a
    // ConZono, whose bounding box is by definition the tightest box enclosing the exact
    // feasible set; since the box-contractor result also encloses that same feasible set,
    // it must contain that tightest bounding box
    auto Z1 = interval_2_zono(b1);
    auto Z2 = interval_2_zono(b2);
    const auto Z_int = intersection(*Z1, *Z2, R);
    const Box zono_bb = Z_int->bounding_box();
    EXPECT_TRUE(result.contains_set(zono_bb));
}

TEST(BoxSetOperations, IntersectionOverDims)
{
    Eigen::Vector<zono_float, -1> lb(3), ub(3);
    lb << 0., 1., 2.; ub << 10., 11., 12.;
    Box b(lb, ub);

    Eigen::Vector<zono_float, -1> lb2(2), ub2(2);
    lb2 << 5., 3.; ub2 << 20., 8.;
    Box other(lb2, ub2);

    const std::vector<int> dims{0, 2};
    Box result = b.intersection_over_dims(other, dims);
    EXPECT_EQ(result.size(), b.size());

    // touched dims are intersected with the corresponding element of other
    EXPECT_TRUE(result.get_element(0) == Interval(5., 10.));
    EXPECT_TRUE(result.get_element(2) == Interval(3., 8.));
    // untouched dim is unaffected
    EXPECT_TRUE(result.get_element(1) == b.get_element(1));

    // definitional: every vertex of the result is in b, and its selected components
    // are in other
    for (const auto& v : box_vertices(result))
    {
        EXPECT_TRUE(b.contains(v));
        Eigen::Vector<zono_float, -1> selected(2);
        selected << v(0), v(2);
        EXPECT_TRUE(other.contains(selected));
    }

    EXPECT_TRUE(intersection_over_dims(b, other, dims) == result);
}

TEST(BoxSetOperations, IntersectionOverDimsDuplicateDims)
{
    Box b(Eigen::Vector<zono_float, -1>::Zero(3), Eigen::Vector<zono_float, -1>::Constant(3, 10.));

    Eigen::Vector<zono_float, -1> lb2(2), ub2(2);
    lb2 << 2., 6.; ub2 << 8., 12.;
    Box other(lb2, ub2);

    // repeated index: dim 1 is intersected cumulatively with both entries of other
    Box result = b.intersection_over_dims(other, {1, 1});
    EXPECT_TRUE(result.get_element(1) == Interval(6., 8.));
    EXPECT_TRUE(result.get_element(0) == b.get_element(0));
    EXPECT_TRUE(result.get_element(2) == b.get_element(2));
}

TEST(BoxSetOperations, IntersectionOverDimsThrows)
{
    Box b(Eigen::Vector<zono_float, -1>::Zero(3), Eigen::Vector<zono_float, -1>::Ones(3));
    Box other(Eigen::Vector<zono_float, -1>::Zero(2), Eigen::Vector<zono_float, -1>::Ones(2));

    EXPECT_THROW(b.intersection_over_dims(other, {-1, 0}), std::invalid_argument);
    EXPECT_THROW(b.intersection_over_dims(other, {3, 0}), std::invalid_argument);
    EXPECT_THROW(b.intersection_over_dims(other, {0}), std::invalid_argument); // dims.size() != other.size()
}
