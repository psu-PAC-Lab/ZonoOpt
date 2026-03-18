#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"
#include <cmath>

using namespace ZonoOpt;

bool check_matrix_equal(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
{
    if (A.rows() != B.rows() || A.cols() != B.cols()) return false;

    for (int i=0; i<A.rows(); ++i) {
        for (int j=0; j<A.cols(); ++j) {
            if (std::abs(A(i,j) - B(i,j)) > zono_eps)
                return false;
        }
    }

    return true;
}

bool check_vector_equal(const Eigen::VectorXd& A, const Eigen::VectorXd& B) {
    if (A.size() != B.size()) return false;

    for (int i=0; i<A.size(); ++i) {
        if (std::abs(A(i,0) - B(i,0)) > zono_eps)
            return false;
    }

    return true;
}

bool check_equal(const HybZono& Z1, HybZono& Z2)
{
    if (Z2.is_0_1_form() != Z1.is_0_1_form())
    {
        Z2.convert_form();
    }

    // make sure dimensions are the same
    if (Z1.get_n() != Z2.get_n() || Z1.get_nGc() != Z2.get_nGc() || Z1.get_nGb() != Z2.get_nGb() || Z1.get_nC() != Z2.get_nC())
        return false;

    // make sure matrices and vectors are the same
    const bool G_equal = check_matrix_equal(Z1.get_G().toDense(), Z2.get_G().toDense());
    const bool c_equal = check_matrix_equal(Z1.get_c(), Z2.get_c());
    const bool A_equal = check_matrix_equal(Z1.get_A(), Z2.get_A());
    const bool b_equal = check_matrix_equal(Z1.get_b(), Z2.get_b());

    return G_equal && c_equal && A_equal && b_equal;
}


int main()
{
    // make 2 zonotopes and a point
    Eigen::VectorXd c1(2);
    c1 << 1., 2.;
    auto Z1 = make_regular_zono_2D(1., 8, false, c1);

    Eigen::VectorXd c2(2);
    c2 << -3., 1.;
    auto Z2 = make_regular_zono_2D(0.1, 6, false, c2);

    Eigen::VectorXd c3(2);
    c3 << 3., -2.;
    Point P3(c3);

    // matrix
    Eigen::MatrixXd M(1,2);
    M << 32., 1.2;
    Eigen::SparseMatrix<double> M_sp = M.sparseView();

    Eigen::MatrixXd M_upper(1, 2);
    M_upper << 33., 1.4;

    IntervalMatrix M_int(M, M_upper);

    // box
    Eigen::VectorXd b1(2);
    Eigen::VectorXd b2(2);
    b1 << 0., 1.;
    b2 << 0.1, 1.04;
    Box box (b1, b2);

    // check operators are consistent with set operations

    // minkowski sum
    auto Z_set = minkowski_sum(*Z1, *Z2);
    auto Z_op = *Z1 + *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum failed");

    // +=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = minkowski_sum(*Z_set, *Z2);
    *Z_op += *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum += failed");

    // minkowski sum with point
    Z_set = minkowski_sum(*Z1, P3);
    Z_op = *Z1 + c3;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with point failed");

    // left sum
    Z_set = minkowski_sum(P3, *Z1);
    Z_op = c3 + *Z1;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with point on left failed");

    // +=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = minkowski_sum(*Z_set, P3);
    *Z_op += c3;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with point += failed");

    // minkowski sum with box
    Z_set = minkowski_sum(*Z1, *interval_2_zono(box));
    Z_op = *Z1 + box;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with box failed");

    // left sum
    Z_set = minkowski_sum(*interval_2_zono(box), *Z1);
    Z_op = box + *Z1;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with box on left failed");

    // +=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = minkowski_sum(*Z_set, *interval_2_zono(box));
    *Z_op += box;
    test_assert(check_equal(*Z_set, *Z_op), "Minkowski sum with box += failed");

    // pontry diff
    Z_set = pontry_diff(*Z1, *Z2, true);
    Z_op = *Z1 - *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference failed");

    // -=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = pontry_diff(*Z_set, *Z2, true);
    *Z_op -= *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference -= failed");

    // pontry diff with point
    Z_set = pontry_diff(*Z1, P3, true);
    Z_op = *Z1 - c3;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference with point failed");

    // -=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = pontry_diff(*Z_set, P3, true);
    *Z_op -= c3;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference with point -= failed");

    // pontry diff with box
    Z_set = pontry_diff(*Z1, *interval_2_zono(box), true);
    Z_op = *Z1 - box;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference with box failed");

    // -=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = pontry_diff(*Z_set, *interval_2_zono(box), true);
    *Z_op -= box;
    test_assert(check_equal(*Z_set, *Z_op), "Pontryagin difference with box -= failed");

    // affine map - sparse
    Z_set = affine_map(*Z1, M_sp);
    Z_op = M_sp * (*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Affine map with sparse matrix failed");

    // affine map - dense
    Z_op = M * (*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Affine map with dense matrix failed");

    // affine inclusion
    Z_set = affine_inclusion(*Z1, M_int);
    Z_op = M_int * (*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Affine inclusion failed");

    // scalar multiplication - left
    double f = 3.2;
    Eigen::SparseMatrix<double> fI = (f*Eigen::MatrixXd::Identity(2, 2)).sparseView();
    Z_set = affine_map(*Z1, fI);
    Z_op = f * (*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Scalar multiplication with left scalar failed");

    // scalar multiplication - right
    Z_op = (*Z1) * f;
    test_assert(check_equal(*Z_set, *Z_op), "Scalar multiplication with right scalar failed");

    // *=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = affine_map(*Z_set, fI);
    *Z_op *= f;
    test_assert(check_equal(*Z_set, *Z_op), "Scalar multiplication *= failed");

    // cartesian product
    Z_set = cartesian_product(*Z1, *Z2);
    Z_op = *Z1 * (*Z2);
    test_assert(check_equal(*Z_set, *Z_op), "Cartesian product failed");

    // *=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = cartesian_product(*Z_set, *Z2);
    *Z_op *= *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Cartesian product *= failed");

    // cartesian product with box
    Z_set = cartesian_product(*Z1, *interval_2_zono(box));
    Z_op = *Z1 * box;
    test_assert(check_equal(*Z_set, *Z_op), "Cartesian product with box failed");

    // cartesian product with box on left
    Z_set = cartesian_product(*interval_2_zono(box), *Z1);
    Z_op = box * (*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Cartesian product with box on left failed");

    // *=
    Z_set.reset(Z1->clone());
    Z_op.reset(Z1->clone());
    Z_set = cartesian_product(*Z_set, *interval_2_zono(box));
    *Z_op *= box;
    test_assert(check_equal(*Z_set, *Z_op), "Cartesian product with box *= failed");

    // intersection
    Z_set = intersection(*Z1, *Z2);
    Z_op = *Z1 & *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Intersection failed");

    // union
    std::vector<std::shared_ptr<HybZono>> Zs;
    Zs.push_back(std::make_shared<HybZono>(*Z1));
    Zs.push_back(std::make_shared<HybZono>(*Z2));
    Z_set = union_of_many(Zs, false, false);
    Z_op = *Z1 | *Z2;
    test_assert(check_equal(*Z_set, *Z_op), "Union failed");

    // unary minus
    Eigen::SparseMatrix<double> mI = (-Eigen::MatrixXd::Identity(2, 2)).sparseView();
    Z_set = affine_map(*Z1, mI);
    Z_op = -(*Z1);
    test_assert(check_equal(*Z_set, *Z_op), "Unary minus failed");

    return 0;
}