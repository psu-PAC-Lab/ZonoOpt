#include <stdexcept>
#include <algorithm>
#include <utility>
#include <memory>
#include <cmath>

#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    std::unique_ptr<HybZono> affine_map(const HybZono& Z,
                                        const Eigen::SparseMatrix<zono_float>& R,
                                        const Eigen::Vector<zono_float, -1>& s)
    {
        // check dimensions
        Eigen::Vector<zono_float, -1> s_def;
        const Eigen::Vector<zono_float, -1>* s_ptr = nullptr;
        if (s.size() == 0) // default argument
        {
            s_def.resize(R.rows());
            s_def.setZero();
            s_ptr = &s_def;
        }
        else
        {
            s_ptr = &s;
        }

        if (R.cols() != Z.n || R.rows() != s_ptr->size())
        {
            throw std::invalid_argument("Linear_map: invalid input dimensions.");
        }

        // early exit
        if (Z.is_empty_set())
        {
            return std::make_unique<EmptySet>(static_cast<int>(R.rows()));
        }

        // apply affine map
        Eigen::SparseMatrix<zono_float> Gc = R * Z.Gc;
        Eigen::SparseMatrix<zono_float> Gb = R * Z.Gb;
        Eigen::Vector<zono_float, -1> c = R * Z.c + *s_ptr;

        // output correct type
        if (Z.is_hybzono())
            return std::make_unique<HybZono>(Gc, Gb, c, Z.Ac, Z.Ab, Z.b, Z.zero_one_form);
        else if (Z.is_conzono())
            return std::make_unique<ConZono>(Gc, c, Z.A, Z.b, Z.zero_one_form);
        else if (Z.is_zono())
            return std::make_unique<Zono>(Gc, c, Z.zero_one_form);
        else
            return std::make_unique<Point>(c);
    }

    std::unique_ptr<HybZono> affine_inclusion(const HybZono& Z, const IntervalMatrix& R,
                                              const Eigen::Vector<zono_float, -1>& s)
    {
        // check dimensions
        Eigen::Vector<zono_float, -1> s_def;
        const Eigen::Vector<zono_float, -1>* s_ptr = nullptr;
        if (s.size() == 0) // default argument
        {
            s_def.resize(R.rows());
            s_def.setZero();
            s_ptr = &s_def;
        }
        else
        {
            s_ptr = &s;
        }

        if (R.cols() != static_cast<size_t>(Z.n) || static_cast<size_t>(s_ptr->size()) != R.rows())
        {
            throw std::invalid_argument("Affine inclusion: invalid input dimensions.");
        }

        // early exit
        if (Z.is_empty_set())
        {
            return std::make_unique<EmptySet>(static_cast<int>(R.rows()));
        }

        // get bounding zonotope
        const auto X_bar = Z.convex_relaxation()->to_zono_approx();
        if (X_bar->is_0_1_form())
        {
            X_bar->convert_form();
        }

        // get m box
        const Box m = R.radius() * X_bar->get_c();

        // get P diagonal matrix
        const Eigen::SparseMatrix<zono_float, Eigen::RowMajor> diam_R_abs_M = R.diam() * X_bar->get_G().cwiseAbs();
        std::vector<Eigen::Triplet<zono_float>> triplets;
        for (size_t i = 0; i < R.rows(); ++i)
        {
            triplets.emplace_back(i, i, 0.5 * (m[i].width() + diam_R_abs_M.row(i).sum()));
        }
        Eigen::SparseMatrix<zono_float> P(R.rows(), R.rows());
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        P.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        P.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // P zonotope
        Zono P_zono(P, *s_ptr, false);

        // affine map
        auto Z_out = affine_map(Z, R.center());

        // minkowski sum P zonotope
        Z_out = minkowski_sum(*Z_out, P_zono);

        return Z_out;
    }

    std::unique_ptr<HybZono> project_onto_dims(const HybZono& Z, const std::vector<int>& dims)
    {
        // make sure all dims are >= 0 and < n
        for (const int dim : dims)
        {
            if (dim < 0 || dim >= Z.n)
            {
                throw std::invalid_argument("Project onto dims: invalid dimension.");
            }
        }

        // early exit
        if (Z.is_empty_set())
        {
            return std::make_unique<EmptySet>(static_cast<int>(dims.size()));
        }

        // build affine map matrix
        Eigen::SparseMatrix<zono_float> R(static_cast<Eigen::Index>(dims.size()), static_cast<Eigen::Index>(Z.n));
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        for (int i = 0; i < static_cast<int>(dims.size()); i++)
        {
            tripvec.emplace_back(i, dims[i], one);
        }
        R.setFromTriplets(tripvec.begin(), tripvec.end());

        // apply affine map
        return affine_map(Z, R);
    }

    std::unique_ptr<HybZono> minkowski_sum(const HybZono& Z1, HybZono& Z2)
    {
        // check dimensions
        if (Z1.n != Z2.n)
        {
            throw std::invalid_argument("Minkowski sum: n dimensions must match.");
        }

        // early exit
        if (Z1.is_empty_set() && Z2.is_empty_set())
        {
            return std::make_unique<EmptySet>(Z1.n);
        }
        else if (Z1.is_empty_set())
        {
            return std::unique_ptr<HybZono>(Z2.clone());
        }
        else if (Z2.is_empty_set())
        {
            return std::unique_ptr<HybZono>(Z1.clone());
        }

        // make sure Z1 and Z2 both using same generator range
        if (Z1.zero_one_form != Z2.zero_one_form)
        {
            Z2.convert_form();
        }

        std::vector<Eigen::Triplet<zono_float>> tripvec;

        Eigen::SparseMatrix<zono_float> Gc = hcat<zono_float>(Z1.Gc, Z2.Gc);
        Eigen::SparseMatrix<zono_float> Gb = hcat<zono_float>(Z1.Gb, Z2.Gb);
        Eigen::Vector<zono_float, -1> c = Z1.c + Z2.c;

        Eigen::SparseMatrix<zono_float> Ac(Z1.nC + Z2.nC, Z1.nGc + Z2.nGc);
        get_triplets_offset<zono_float>(Z1.Ac, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ac, tripvec, Z1.nC, Z1.nGc);
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ab(Z1.nC + Z2.nC, Z1.nGb + Z2.nGb);
        get_triplets_offset<zono_float>(Z1.Ab, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ab, tripvec, Z1.nC, Z1.nGb);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());

        Eigen::Vector<zono_float, -1> b(Z1.nC + Z2.nC);
        b.segment(0, Z1.nC) = Z1.b;
        b.segment(Z1.nC, Z2.nC) = Z2.b;

        // return correct output type
        if (Z1.is_hybzono() || Z2.is_hybzono())
            return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, Z1.zero_one_form);
        else if (Z1.is_conzono() || Z2.is_conzono())
            return std::make_unique<ConZono>(Gc, c, Ac, b, Z1.zero_one_form);
        else
            return std::make_unique<Zono>(Gc, c, Z1.zero_one_form);
    }

    std::unique_ptr<HybZono> intersection(const HybZono& Z1, HybZono& Z2, const Eigen::SparseMatrix<zono_float>& R)
    {
        // handle default arguments
        const Eigen::SparseMatrix<zono_float>* R_ptr = nullptr;
        Eigen::SparseMatrix<zono_float> R_def;
        if (R.rows() == 0 && R.cols() == 0)
        {
            R_def.resize(Z1.n, Z1.n);
            R_def.setIdentity();
            R_ptr = &R_def;
        }
        else
        {
            R_ptr = &R;
        }

        // check dimensions
        if (R_ptr->rows() != Z2.n || R_ptr->cols() != Z1.n)
        {
            throw std::invalid_argument("Intersection: inconsistent input dimensions.");
        }

        // early exit
        if (Z1.is_empty_set() || Z2.is_empty_set())
        {
            return std::make_unique<EmptySet>(Z1.n);
        }

        // make sure Z1 and Z2 both using same generator range
        if (Z1.zero_one_form != Z2.zero_one_form)
        {
            Z2.convert_form();
        }

        // compute intersection
        Eigen::SparseMatrix<zono_float> Gc = Z1.Gc;
        Gc.conservativeResize(Z1.n, Z1.nGc + Z2.nGc);

        Eigen::SparseMatrix<zono_float> Gb = Z1.Gb;
        Gb.conservativeResize(Z1.n, Z1.nGb + Z2.nGb);

        Eigen::Vector<zono_float, -1> c = Z1.c;

        std::vector<Eigen::Triplet<zono_float>> tripvec;
        Eigen::SparseMatrix<zono_float> Ac(Z1.nC + Z2.nC + R_ptr->rows(), Z1.nGc + Z2.nGc);
        get_triplets_offset<zono_float>(Z1.Ac, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ac, tripvec, Z1.nC, Z1.nGc);
        Eigen::SparseMatrix<zono_float> RZ1Gc = (*R_ptr) * Z1.Gc;
        get_triplets_offset<zono_float>(RZ1Gc, tripvec, Z1.nC + Z2.nC, 0);
        Eigen::SparseMatrix<zono_float> mZ2Gc = -Z2.Gc;
        get_triplets_offset<zono_float>(mZ2Gc, tripvec, Z1.nC + Z2.nC, Z1.nGc);
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ab(Z1.nC + Z2.nC + R_ptr->rows(), Z1.nGb + Z2.nGb);
        get_triplets_offset<zono_float>(Z1.Ab, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ab, tripvec, Z1.nC, Z1.nGb);
        Eigen::SparseMatrix<zono_float> RZ1Gb = (*R_ptr) * Z1.Gb;
        get_triplets_offset<zono_float>(RZ1Gb, tripvec, Z1.nC + Z2.nC, 0);
        Eigen::SparseMatrix<zono_float> mZ2Gb = -Z2.Gb;
        get_triplets_offset<zono_float>(mZ2Gb, tripvec, Z1.nC + Z2.nC, Z1.nGb);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());


        Eigen::Vector<zono_float, -1> b(Z1.nC + Z2.nC + R_ptr->rows());
        b.segment(0, Z1.nC) = Z1.b;
        b.segment(Z1.nC, Z2.nC) = Z2.b;
        b.segment(Z1.nC + Z2.nC, R_ptr->rows()) = Z2.c - (*R_ptr) * Z1.c;

        // return correct output type
        if (Z1.is_hybzono() || Z2.is_hybzono())
            return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, Z1.zero_one_form);
        else
            return std::make_unique<ConZono>(Gc, c, Ac, b, Z1.zero_one_form);
    }

    std::unique_ptr<HybZono> intersection_over_dims(const HybZono& Z1,
                                                    HybZono& Z2, const std::vector<int>& dims)
    {
        // check dimensions
        if (static_cast<size_t>(Z2.n) != dims.size())
        {
            throw std::invalid_argument("Intersection over dims: Z2.n must equal number of dimensions.");
        }

        // make sure dims are >=0 and <Z1.n
        for (const int dim : dims)
        {
            if (dim < 0 || dim >= Z1.n)
            {
                throw std::invalid_argument("Intersection over dims: invalid dimension.");
            }
        }

        // build projection matrix
        Eigen::SparseMatrix<zono_float> R(static_cast<Eigen::Index>(dims.size()), static_cast<Eigen::Index>(Z1.n));
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        for (int i = 0; i < static_cast<int>(dims.size()); ++i)
        {
            tripvec.emplace_back(i, dims[i], one);
        }
        R.setFromTriplets(tripvec.begin(), tripvec.end());

        // generalized intersection
        return intersection(Z1, Z2, R);
    }

    std::unique_ptr<HybZono> halfspace_intersection(HybZono& Z, const Eigen::SparseMatrix<zono_float>& H,
                                                    const Eigen::Vector<zono_float, -1>& f,
                                                    const Eigen::SparseMatrix<zono_float>& R)
    {
        // use the constrain function

        // convert H to row-major to efficiently convert into Inequality form
        Eigen::SparseMatrix<zono_float, Eigen::RowMajor> H_rm = H;
        std::vector<Inequality> ineqs;
        for (int k = 0; k < H_rm.rows(); ++k)
        {
            // get row of constraint matrix
            std::vector<Eigen::Triplet<zono_float>> trips_row = get_triplets_row<zono_float>(H_rm, k);

            // build constraint
            Inequality ineq(static_cast<int>(H_rm.cols()));
            for (const auto& trip : trips_row)
            {
                ineq.add_term(trip.col(), trip.value());
            }
            ineq.set_rhs(f(k));
            ineq.set_ineq_type(LESS_OR_EQUAL);

            ineqs.push_back(ineq);
        }

        // call constrain
        return constrain(Z, ineqs, R);
    }

    // pontryagin difference
    std::unique_ptr<HybZono> pontry_diff(HybZono& Z1, HybZono& Z2, bool exact)
    {
        // check dimensions
        if (Z1.n != Z2.n)
        {
            throw std::invalid_argument("Pontryagin difference: n dimensions must match.");
        }

        // check empty set inputs
        if (Z1.is_empty_set() || Z2.is_empty_set())
            return std::unique_ptr<HybZono>(Z1.clone());

        // check point inputs
        if (Z1.nG == 0 && !exact)
            return std::make_unique<EmptySet>(Z1.n);
        if (Z2.nG == 0)
            exact = true; // easy to compute exactly

        // check no exact pontry diff with conzono subtrahend
        if (exact && (Z2.is_conzono() || Z2.is_hybzono()))
            throw std::invalid_argument(
                "Pontryagin difference: cannot compute exact set when subtrahend is a ConZono or HybZono");

        // require Z1 and Z2 to be in [-1,1] form
        if (Z1.zero_one_form) Z1.convert_form();
        if (Z2.zero_one_form) Z2.convert_form();

        // exact case
        if (exact)
        {
            // init Zout
            std::unique_ptr<HybZono> Z_out(Z1.clone());
            Z_out->c -= Z2.c;

            // iteratively compute pontryagin difference from columns of Z2 generator matrix
            const Eigen::Matrix<zono_float, -1, -1> G2 = Z2.G.toDense();

            for (int i = 0; i < Z2.nG; ++i)
            {
                std::unique_ptr<HybZono> Z_plus(Z_out->clone()), Z_minus(Z_out->clone());
                Z_plus->c += G2.col(i);
                Z_minus->c -= G2.col(i);
                Z_out = intersection(*Z_plus, *Z_minus);
            }

            return Z_out;
        }

        // inexact case

        // logic to handle hybrid zonotope cases. Avoiding recursion to prevent redundant get_leaves calculations.
        if (Z1.is_hybzono() && Z2.is_hybzono())
        {
            // get leaves and convert to zonos
            auto Z1_leaves = Z1.get_leaves();
            auto Z2_CZ_leaves = Z2.get_leaves();
            std::vector<Zono> Z2_leaves;
            for (auto& CZ : Z2_CZ_leaves)
            {
                Z2_leaves.push_back(*CZ.to_zono_approx());
            }

            // take union of pontry diffs for inner approx
            std::vector<std::shared_ptr<HybZono>> leaf_diffs; // init
            for (auto& Z1_leaf : Z1_leaves)
            {
                std::unique_ptr<HybZono> Z_out; // declare
                for (auto& Z2_leaf : Z2_leaves)
                {
                    if (!Z_out)
                    {
                        Z_out = pontry_diff(Z1_leaf, Z2_leaf, false);
                    }
                    else
                    {
                        Z_out = intersection(*Z_out, *pontry_diff(Z1_leaf, Z2_leaf, false));
                    }
                }
                leaf_diffs.push_back(std::move(Z_out));
            }
            return union_of_many(leaf_diffs);
        }
        else if (Z1.is_hybzono())
        {
            // get leaves and convert to zonos
            auto Z1_leaves = Z1.get_leaves();
            auto Z2_CZ = dynamic_cast<const ConZono*>(&Z2);
            Zono Z2_zono = *Z2_CZ->to_zono_approx();

            // take union of pontry diffs for inner approx
            std::vector<std::shared_ptr<HybZono>> leaf_diffs; // init
            for (auto& Z1_leaf : Z1_leaves)
            {
                leaf_diffs.push_back(pontry_diff(Z1_leaf, Z2_zono, false));
            }
            return union_of_many(leaf_diffs);
        }
        else if (Z2.is_hybzono())
        {
            // get leaves and convert to zonos
            auto Z2_CZ_leaves = Z2.get_leaves();
            std::vector<Zono> Z2_leaves;
            for (auto& CZ : Z2_CZ_leaves)
            {
                Z2_leaves.push_back(*CZ.to_zono_approx());
            }

            // take intersection of pontry diffs
            std::unique_ptr<HybZono> Z_out; // declare
            for (auto& Z2_leaf : Z2_leaves)
            {
                if (!Z_out)
                {
                    Z_out = pontry_diff(Z1, Z2_leaf, false);
                }
                else
                {
                    Z_out = intersection(*Z_out, *pontry_diff(Z1, Z2_leaf, false));
                }
            }
            return Z_out;
        }
        else
        {
            // convert to zono subtrahend (do nothing if already zono)
            auto Z2_CZ = dynamic_cast<const ConZono*>(&Z2);
            Zono Z2_zono = *Z2_CZ->to_zono_approx();

            // get [G; A] matrices
            const Eigen::SparseMatrix<zono_float> GA1 = vcat<zono_float>(Z1.G, Z1.A);
            Eigen::SparseMatrix<zono_float> GA2 = Z2_zono.G;
            GA2.conservativeResize(Z2_zono.n + Z1.nC, Z2_zono.nG);

            Eigen::SparseMatrix<zono_float> M;
            if (Z1.n + Z1.nC == Z1.nG)
            {
                // solve GA1^{-1}*GA2
                Eigen::SparseLU<Eigen::SparseMatrix<zono_float>> lu(GA1);
                if (lu.info() != Eigen::Success)
                {
                    throw std::runtime_error(
                        "Pontryagin difference: failed to peform LU decomposition. Most likely [G; A] is not full row rank");
                }
                M = lu.solve(GA2);
            }
            else
            {
                // solve pinv(GA1)*GA2 where pinv is right pseudoinverse
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>> ldlt(GA1 * GA1.transpose());
                if (ldlt.info() != Eigen::Success)
                {
                    throw std::runtime_error(
                        "Pontryagin difference: failed to perform LDLT decomposition. Most likely [G; A] is not full row rank");
                }
                const Eigen::SparseMatrix<zono_float> ldlt_sol = ldlt.solve(GA2);
                M = GA1.transpose() * ldlt_sol;
            }

            std::vector<Eigen::Triplet<zono_float>> triplets;
            triplets.reserve(Z1.nG);
            Eigen::Vector<zono_float, -1> e(Z1.nG);
            for (int i = 0; i < Z1.nG; ++i)
            {
                e.setZero();
                e(i) = 1.0;
                const zono_float d = 1 - (e.transpose() * M).cwiseAbs().sum();
                if (d < 0)
                {
                    return std::make_unique<EmptySet>(Z1.n);
                }
                triplets.emplace_back(i, i, d);
            }
            Eigen::SparseMatrix<zono_float> D(Z1.nG, Z1.nG);
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
            D.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
            D.setFromTriplets(triplets.begin(), triplets.end());
#endif

            return std::make_unique<ConZono>(Z1.G * D, Z1.c - Z2_zono.c, Z1.A * D, Z1.b, false);
        }
    }

    std::unique_ptr<HybZono> union_of_many(const std::vector<std::shared_ptr<HybZono>>& Zs_in,
                                           const bool preserve_sharpness, const bool expose_indicators)
    {
        // remove empty sets from sets to be unioned
        std::vector<std::shared_ptr<HybZono>> Zs;
        for (auto& Z : Zs_in)
        {
            if (!Z->is_empty_set())
            {
                Zs.push_back(Z);
            }
        }

        // check we are taking a union of at least one zonotope
        if (Zs.empty())
        {
            throw std::invalid_argument("Union: empty input vector.");
        }

        // check dimensions
        const int n = Zs[0]->n;
        for (const auto& Z : Zs)
        {
            if (Z->n != n)
            {
                throw std::invalid_argument("Union: inconsistent dimensions.");
            }
        }

        // make sure all Zs are using [0,1] generators
        for (const auto& Z : Zs)
        {
            if (!Z->zero_one_form)
            {
                Z->convert_form();
            }
        }

        // declare
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        int rows = 0, cols = 0;
        std::vector<int> idx_sum_to_1;
        Eigen::SparseMatrix<zono_float> Gc, Gb, Ac, Ab;
        Eigen::Vector<zono_float, -1> c, b;

        if (preserve_sharpness)
        {
            // constraints

            // Ac
            tripvec.clear();
            rows = 0;
            cols = 0;
            for (const auto& Z : Zs)
            {
                // equality constraints
                get_triplets_offset<zono_float>(Z->Ac, tripvec, rows, cols);
                rows += Z->nC;

                // identity matrices
                for (int i = 0; i < Z->nGc; i++)
                {
                    tripvec.emplace_back(rows + i, cols + i, one);
                }
                for (int i = 0; i < Z->nG; i++)
                {
                    tripvec.emplace_back(rows + i, cols + Z->nGc + i, one);
                }

                // increment
                rows += Z->nG;
                cols += Z->nGc + Z->nG;
            }
            Ac.resize(rows + 1, cols); // last row all zeroes
            Ac.setFromTriplets(tripvec.begin(), tripvec.end());

            // Ab
            tripvec.clear();
            rows = 0;
            cols = 0;
            for (const auto& Z : Zs)
            {
                // equality constraints
                get_triplets_offset<zono_float>(Z->Ab, tripvec, rows, cols);

                // last column
                for (int i = 0; i < Z->nC; i++)
                {
                    tripvec.emplace_back(rows + i, cols + Z->nGb, -Z->b(i));
                }
                rows += Z->nC;

                // identity matrix
                for (int i = 0; i < Z->nGb; i++)
                {
                    tripvec.emplace_back(rows + Z->nGc + i, cols + i, one);
                }

                // last column
                for (int i = 0; i < Z->nG; i++)
                {
                    tripvec.emplace_back(rows + i, cols + Z->nGb, -one);
                }

                // increment
                rows += Z->nG;
                cols += Z->nGb + 1;

                // track sum-to-1 binaries
                idx_sum_to_1.push_back(cols - 1);
            }
            // sum to 1 constraint
            for (int& it : idx_sum_to_1)
            {
                tripvec.emplace_back(rows, it, one);
            }
            rows++;
            Ab.resize(rows, cols);
            Ab.setFromTriplets(tripvec.begin(), tripvec.end());

            // b
            b.resize(Ab.rows());
            b.setZero();
            b(Ab.rows() - 1) = 1.0;

            // generators
            int n_out = n;
            if (expose_indicators)
                n_out += static_cast<int>(Zs.size());

            // Gc
            tripvec.clear();
            cols = 0;
            for (const auto& Z : Zs)
            {
                get_triplets_offset<zono_float>(Z->Gc, tripvec, 0, cols);
                cols += 2 * Z->nGc + Z->nGb;
            }
            Gc.resize(n_out, cols);
            Gc.setFromTriplets(tripvec.begin(), tripvec.end());

            // Gb
            tripvec.clear();
            cols = 0;
            for (const auto& Z : Zs)
            {
                get_triplets_offset<zono_float>(Z->Gb, tripvec, 0, cols);
                cols += Z->nGb;
                for (int i = 0; i < Z->n; i++)
                {
                    tripvec.emplace_back(i, cols, Z->c(i));
                }
                cols += 1;
            }
            for (int i = n; i < n_out; ++i)
            {
                tripvec.emplace_back(i, idx_sum_to_1[i - n], one);
            }
            Gb.resize(n_out, cols);
            Gb.setFromTriplets(tripvec.begin(), tripvec.end());

            // c
            c.resize(n_out);
            c.setZero();
        }
        else
        {
            // constraints

            // Ac
            tripvec.clear();
            rows = 0;
            cols = 0;
            for (const auto& Z : Zs)
            {
                // populate first row
                for (int i = 0; i < Z->nGc; i++)
                {
                    tripvec.emplace_back(rows, cols + i, one);
                }
                tripvec.emplace_back(rows, cols + Z->nGc, static_cast<zono_float>(Z->nG));

                // increment
                rows++;

                // equality constraints
                get_triplets_offset<zono_float>(Z->Ac, tripvec, rows, cols);

                // increment
                rows += Z->nC;
                cols += Z->nGc + 1;
            }
            Ac.resize(rows + 1, cols); // last row all zeroes
            Ac.setFromTriplets(tripvec.begin(), tripvec.end());

            // Ab
            tripvec.clear();
            rows = 0;
            cols = 0;
            for (const auto& Z : Zs)
            {
                // populate first row
                for (int i = 0; i < Z->nGb; i++)
                {
                    tripvec.emplace_back(rows, cols + i, one);
                }
                tripvec.emplace_back(rows, cols + Z->nGb, static_cast<zono_float>(-Z->nG));

                // increment
                rows++;

                // equality constraints
                get_triplets_offset<zono_float>(Z->Ab, tripvec, rows, cols);

                // last column
                for (int i = 0; i < Z->nC; i++)
                {
                    tripvec.emplace_back(rows + i, cols + Z->nGb, -Z->b(i));
                }

                // increment
                rows += Z->nC;
                cols += Z->nGb + 1;

                // track sum-to-1 binaries
                idx_sum_to_1.push_back(cols - 1);
            }
            // sum to 1 constraint
            for (int& it : idx_sum_to_1)
            {
                tripvec.emplace_back(rows, it, one);
            }
            rows++;
            Ab.resize(rows, cols);
            Ab.setFromTriplets(tripvec.begin(), tripvec.end());

            // b
            b.resize(Ab.rows());
            b.setZero();
            b(Ab.rows() - 1) = 1.0;

            // generators
            int n_out = n;
            if (expose_indicators)
                n_out += static_cast<int>(Zs.size());

            // Gc
            tripvec.clear();
            cols = 0;
            for (const auto& Z : Zs)
            {
                get_triplets_offset<zono_float>(Z->Gc, tripvec, 0, cols);
                cols += Z->nGc + 1;
            }
            Gc.resize(n_out, cols);
            Gc.setFromTriplets(tripvec.begin(), tripvec.end());

            // Gb
            tripvec.clear();
            cols = 0;
            for (const auto& Z : Zs)
            {
                get_triplets_offset<zono_float>(Z->Gb, tripvec, 0, cols);
                cols += Z->nGb;
                for (int i = 0; i < Z->n; i++)
                {
                    tripvec.emplace_back(i, cols, Z->c(i));
                }
                cols += 1;
            }
            for (int i = n; i < n_out; ++i)
            {
                tripvec.emplace_back(i, idx_sum_to_1[i - n], one);
            }
            Gb.resize(n_out, cols);
            Gb.setFromTriplets(tripvec.begin(), tripvec.end());

            // c
            c.resize(n_out);
            c.setZero();
        }

        // check if known to be sharp
        bool sharp = preserve_sharpness;
        size_t i = 0;
        while (sharp && i < Zs.size())
        {
            sharp = Zs[i]->sharp;
            ++i;
        }

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, true, sharp);
    }

    std::unique_ptr<ConZono> convex_hull(const std::vector<std::shared_ptr<HybZono>>& Zs_in)
    {
        // make sure all sets are sharp
        for (const auto& Z : Zs_in)
        {
            if (!Z->sharp)
            {
                throw std::invalid_argument("Convex hull: all input sets must be sharp.");
            }
        }

        // get union
        const std::unique_ptr<HybZono> Z_union = union_of_many(Zs_in, true);

        // take convex relaxation
        return Z_union->convex_relaxation();
    }

    std::unique_ptr<HybZono> cartesian_product(const HybZono& Z1, HybZono& Z2)
    {
        // trivial case
        if (Z1.is_empty_set() || Z2.is_empty_set())
        {
            return std::make_unique<EmptySet>(Z1.n + Z2.n);
        }

        // make sure Z1 and Z2 both using same generator range
        if (Z1.zero_one_form != Z2.zero_one_form)
        {
            Z2.convert_form();
        }

        // declare
        std::vector<Eigen::Triplet<zono_float>> tripvec;

        // take Cartesian product
        Eigen::SparseMatrix<zono_float> Gc(Z1.n + Z2.n, Z1.nGc + Z2.nGc);
        get_triplets_offset<zono_float>(Z1.Gc, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Gc, tripvec, Z1.n, Z1.nGc);
        Gc.setFromTriplets(tripvec.begin(), tripvec.end());

        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Gb(Z1.n + Z2.n, Z1.nGb + Z2.nGb);
        get_triplets_offset<zono_float>(Z1.Gb, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Gb, tripvec, Z1.n, Z1.nGb);
        Gb.setFromTriplets(tripvec.begin(), tripvec.end());

        Eigen::Vector<zono_float, -1> c(Z1.n + Z2.n);
        c.segment(0, Z1.n) = Z1.c;
        c.segment(Z1.n, Z2.n) = Z2.c;

        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ac(Z1.nC + Z2.nC, Z1.nGc + Z2.nGc);
        get_triplets_offset<zono_float>(Z1.Ac, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ac, tripvec, Z1.nC, Z1.nGc);
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ab(Z1.nC + Z2.nC, Z1.nGb + Z2.nGb);
        get_triplets_offset<zono_float>(Z1.Ab, tripvec, 0, 0);
        get_triplets_offset<zono_float>(Z2.Ab, tripvec, Z1.nC, Z1.nGb);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());

        Eigen::Vector<zono_float, -1> b(Z1.nC + Z2.nC);
        b.segment(0, Z1.nC) = Z1.b;
        b.segment(Z1.nC, Z2.nC) = Z2.b;

        // return correct output type
        if (Z1.is_hybzono() || Z2.is_hybzono())
            return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, Z1.zero_one_form);
        else if (Z1.is_conzono() || Z2.is_conzono())
            return std::make_unique<ConZono>(Gc, c, Ac, b, Z1.zero_one_form);
        else if (Z1.is_zono() || Z2.is_zono())
            return std::make_unique<Zono>(Gc, c, Z1.zero_one_form);
        else
            return std::make_unique<Point>(c);
    }

    std::unique_ptr<HybZono> constrain(HybZono& Z, const std::vector<Inequality>& ineqs,
                                       const Eigen::SparseMatrix<zono_float>& R)
    {
        // trivial case
        if (Z.is_empty_set())
        {
            return std::make_unique<EmptySet>(Z.n);
        }

        // handle default arguments
        const Eigen::SparseMatrix<zono_float>* R_ptr = nullptr;
        Eigen::SparseMatrix<zono_float> R_def;
        if (R.rows() == 0 && R.cols() == 0)
        {
            R_def.resize(Z.n, Z.n);
            R_def.setIdentity();
            R_ptr = &R_def;
        }
        else
        {
            R_ptr = &R;
        }

        // check that dimensions match
        for (const auto& ineq : ineqs)
        {
            if (R_ptr->rows() != ineq.get_n_dims())
                throw std::invalid_argument("Inequality does not have the same number of dimensions as set");
        }

        // make sure Z in 0-1 form
        if (!Z.is_0_1_form())
            Z.convert_form();

        // build constraints
        std::vector<Eigen::Triplet<zono_float>> triplets_Ac, triplets_Ab;
        Eigen::Vector<zono_float, -1> b_new(ineqs.size());
        Eigen::Index n_cons = 0, n_slack = 0;
        Eigen::SparseMatrix<zono_float, Eigen::RowMajor> RGc = (*R_ptr) * Z.Gc;
        Eigen::SparseMatrix<zono_float, Eigen::RowMajor> RGb = (*R_ptr) * Z.Gb;
        Eigen::Vector<zono_float, -1> Rc = (*R_ptr) * Z.c;

        for (const auto& ineq : ineqs)
        {
            zono_float rhs;
            switch (ineq.get_ineq_type())
            {
            case LESS:
                rhs = ineq.get_rhs() - zono_eps;
                break;
            case GREATER:
                rhs = ineq.get_rhs() + zono_eps;
                break;
            default:
                rhs = ineq.get_rhs();
                break;
            }
            zono_float gamma = rhs; // slack variable scaling
            zono_float db = 0;
            const auto ineq_type = ineq.get_ineq_type();

            for (const auto& [idx, coeff] : ineq.get_terms())
            {
                const auto trips_Gc = get_triplets_row<zono_float>(RGc, idx);
                for (const auto& trip : trips_Gc)
                {
                    const zono_float val = coeff * trip.value();
                    triplets_Ac.emplace_back(static_cast<int>(n_cons), static_cast<int>(trip.col()), val);
                    if ((val < 0 && (ineq_type == LESS_OR_EQUAL || ineq_type == LESS)) || (val > 0 && (ineq_type ==
                        GREATER_OR_EQUAL || ineq_type == GREATER)))
                        gamma -= val;
                }

                const auto trips_Gb = get_triplets_row<zono_float>(RGb, idx);
                for (const auto& trip : trips_Gb)
                {
                    const zono_float val = coeff * trip.value();
                    triplets_Ab.emplace_back(static_cast<int>(n_cons), static_cast<int>(trip.col()), val);
                    if ((val < 0 && (ineq_type == LESS_OR_EQUAL || ineq_type == LESS)) || (val > 0 && (ineq_type ==
                        GREATER_OR_EQUAL || ineq_type == GREATER)))
                        gamma -= val;
                }

                const zono_float ddb = coeff * Rc(idx);
                db -= ddb;
                gamma -= ddb;
            }

            // rhs
            b_new(n_cons) = rhs + db;

            // add slack variable
            if (ineq_type != EQUAL)
            {
                triplets_Ac.emplace_back(static_cast<int>(n_cons), Z.nGc + static_cast<int>(n_slack), gamma);
                ++n_slack; // increment slack variable index
            }

            // increment
            ++n_cons;
        }

        Eigen::SparseMatrix<zono_float> Ac_cons(static_cast<Eigen::Index>(ineqs.size()), Z.nGc + n_slack);
        Eigen::SparseMatrix<zono_float> Ab_cons(static_cast<Eigen::Index>(ineqs.size()), Z.nGb);
        Ac_cons.setFromTriplets(triplets_Ac.begin(), triplets_Ac.end());
        Ab_cons.setFromTriplets(triplets_Ab.begin(), triplets_Ab.end());

        // set matrices / vectors
        Eigen::SparseMatrix<zono_float> Z_Ac = Z.Ac;
        Z_Ac.conservativeResize(Z.nC, Z.nGc + n_slack);
        Eigen::SparseMatrix<zono_float> Ac = vcat(Z_Ac, Ac_cons);
        Eigen::SparseMatrix<zono_float> Ab = vcat(Z.Ab, Ab_cons);
        Eigen::Vector<zono_float, -1> b(Z.nC + static_cast<Eigen::Index>(ineqs.size()));
        b.segment(0, Z.nC) = Z.b;
        b.segment(Z.nC, static_cast<Eigen::Index>(ineqs.size())) = b_new;

        Eigen::SparseMatrix<zono_float> Z_Gc = Z.Gc;
        Z_Gc.conservativeResize(Z.n, Z.nGc + n_slack);

        // output correct type
        if (Z.is_hybzono())
            return std::make_unique<HybZono>(Z_Gc, Z.Gb, Z.c, Ac, Ab, b, true);
        else
            return std::make_unique<ConZono>(Z_Gc, Z.c, Ac, b, true);
    }

    std::unique_ptr<HybZono> set_diff(const HybZono& Z1, HybZono& Z2, const zono_float delta_m,
                                      const bool remove_redundancy,
                                      const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                                      const int n_leaves,
                                      const int contractor_iter)
    {
        // trivial case
        if (Z2.is_empty_set())
        {
            return std::unique_ptr<HybZono>(Z1.clone());
        }

        // get complement of Z2
        auto Z2_comp = Z2.complement(delta_m, remove_redundancy, settings, solution, n_leaves, contractor_iter);

        // set difference
        return intersection(Z1, *Z2_comp);
    }

    std::unique_ptr<HybZono> zono_union_2_hybzono(std::vector<std::shared_ptr<Zono>>& Zs, const bool expose_indicators)
    {
        // can't be empty
        if (Zs.empty())
        {
            throw std::invalid_argument("Zono union: empty input vector.");
        }

        // zonotope dimension
        int n_dims = Zs[0]->n;
        const int n_zonos = static_cast<int>(Zs.size());

        // loop through Zs
        for (auto& Z : Zs)
        {
            // make sure dimensions are consistent
            if (Z->n != n_dims)
            {
                throw std::invalid_argument("Zono union: inconsistent dimensions.");
            }

            // convert to [0,1] form
            if (!Z->zero_one_form)
            {
                Z->convert_form();
            }
        }

        // get unique generators and incidence matrix

        // initialize S and M matrices as std::vectors
        // each entry is a row
        int n_gens;
        std::vector<Eigen::Matrix<zono_float, 1, -1>> M_vec;
        std::vector<Eigen::Matrix<zono_float, -1, 1>> S_vec;
        Eigen::Matrix<zono_float, 1, -1> M_row(n_zonos);
        Eigen::Matrix<zono_float, -1, -1> Gd;

        // loop through each polytope
        for (int i = 0; i < n_zonos; i++)
        {
            n_gens = Zs[i]->nG;
            Gd = Zs[i]->G.toDense();
            for (int j = 0; j < n_gens; j++)
            {
                // check if the generator is already in S_vec
                auto generator_equal = [&](const Eigen::Matrix<zono_float, -1, 1>& s) -> bool
                {
                    return (s - Gd.col(j)).norm() < zono_eps;
                };

                if (auto it_S = std::find_if(S_vec.begin(), S_vec.end(), generator_equal); it_S == S_vec.end())
                {
                    S_vec.emplace_back(Gd.col(j));
                    M_row.setZero();
                    M_row(i) = 1;
                    M_vec.push_back(M_row);
                }
                else
                {
                    const int idx = static_cast<int>(std::distance(S_vec.begin(), it_S));
                    M_vec[idx](i) = 1;
                }
            }
        }

        const int nG = static_cast<int>(S_vec.size()); // number of unique generators

        // convert to Eigen matrices
        Eigen::Matrix<zono_float, -1, -1> S(n_dims, nG);
        Eigen::Matrix<zono_float, -1, -1> M(nG, n_zonos);
        for (int i = 0; i < nG; i++)
        {
            S.col(i) = S_vec[i];
            M.row(i) = M_vec[i];
        }

        // directly build hybzono in [0,1] form

        // declare
        std::vector<Eigen::Triplet<zono_float>> tripvec;

        // output dimension
        int n_out = n_dims;
        if (expose_indicators)
            n_out += static_cast<int>(n_zonos);

        // Gc = [S, 0]
        Eigen::SparseMatrix<zono_float> Gc = S.sparseView();
        Gc.conservativeResize(n_out, 2 * nG);

        // Gb = [c0, c1, ...]
        tripvec.clear();
        for (int i = 0; i < n_zonos; ++i)
        {
            for (int j = 0; j < n_dims; ++j)
            {
                tripvec.emplace_back(j, i, Zs[i]->c(j));
            }
        }
        for (int i = n_dims; i < n_out; ++i)
        {
            tripvec.emplace_back(i, i - n_dims, one);
        }
        Eigen::SparseMatrix<zono_float> Gb(n_out, n_zonos);
        Gb.setFromTriplets(tripvec.begin(), tripvec.end());

        // c = 0
        Eigen::Vector<zono_float, -1> c(n_out);
        c.setZero();

        // Ac = [0^T, 0^T;
        //       I, diag[sum(M, 2)]]
        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ac(1 + nG, 2 * nG);
        Eigen::SparseMatrix<zono_float> I_ng(nG, nG);
        I_ng.setIdentity();
        get_triplets_offset<zono_float>(I_ng, tripvec, 1, 0);
        Eigen::Vector<zono_float, -1> sum_M = M.rowwise().sum();
        for (int i = 0; i < nG; i++)
        {
            tripvec.emplace_back(1 + i, nG + i, sum_M(i));
        }
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        // Ab = [1^T;
        //       -M]
        Eigen::SparseMatrix<zono_float> Ab(1 + nG, n_zonos);
        tripvec.clear();
        for (int i = 0; i < n_zonos; i++)
        {
            tripvec.emplace_back(0, i, one);
        }
        Eigen::SparseMatrix<zono_float> mM_sp = -M.sparseView();
        get_triplets_offset<zono_float>(mM_sp, tripvec, 1, 0);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());

        // b = [1;
        //      0]
        Eigen::Vector<zono_float, -1> b(1 + nG);
        b.setZero();
        b(0) = 1;

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, true, true);
    }

}