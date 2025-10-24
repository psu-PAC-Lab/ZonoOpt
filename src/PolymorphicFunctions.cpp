/**
 * @file PolymorphicFunctions.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Set operations and setup functions for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <stdexcept>
#include <algorithm>
#include <utility>
#include <memory>
#include <cmath>

#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    // type checking
    bool HybZono::is_point() const
    {
        const auto PointCast = dynamic_cast<const Point*>(this);
        return PointCast != nullptr;
    }

    bool HybZono::is_zono() const
    {
        const auto ZonoCast = dynamic_cast<const Zono*>(this);
        const auto PointCast = dynamic_cast<const Point*>(this);
        return (ZonoCast != nullptr) && (PointCast == nullptr);
    }

    bool HybZono::is_conzono() const
    {
        const auto ConZonoCast = dynamic_cast<const ConZono*>(this);
        const auto ZonoCast = dynamic_cast<const Zono*>(this);
        const auto EmptySetCast = dynamic_cast<const EmptySet*>(this);

        return (ConZonoCast != nullptr) && (ZonoCast == nullptr) && (EmptySetCast == nullptr);
    }

    bool HybZono::is_hybzono() const
    {
        const auto HybZonoCast = dynamic_cast<const HybZono*>(this);
        const auto ConZonoCast = dynamic_cast<const ConZono*>(this);

        return (HybZonoCast != nullptr) && (ConZonoCast == nullptr);
    }

    bool HybZono::is_empty_set() const
    {
        const auto EmptySetCast = dynamic_cast<const EmptySet*>(this);
        return EmptySetCast != nullptr;
    }


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
#if EIGEN_VERSION_AT_LEAST(5,0,0)
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

    std::unique_ptr<HybZono> HybZono::do_complement(const zono_float delta_m, const bool remove_redundancy,
                                                    const OptSettings& settings,
                                                    OptSolution* solution, const int n_leaves,
                                                    const int contractor_iter)
    {
        // make sure set in [-1,1] form
        if (this->is_0_1_form()) this->convert_form();

        // need to get leaves and do complement for each leaf if Z is a hybzono
        auto leaves = this->get_leaves(remove_redundancy, settings, solution, n_leaves, contractor_iter);
        if (leaves.empty())
        {
            throw std::runtime_error("HybZono complement: set is empty.");
        }
        std::vector<std::unique_ptr<HybZono>> complements; // init
        for (auto& leaf : leaves)
        {
            complements.emplace_back(leaf.complement(delta_m));
        }
        std::unique_ptr<HybZono> Z_out;
        for (auto& comp : complements)
        {
            if (!Z_out)
            {
                Z_out.reset(comp->clone());
            }
            else
            {
                Z_out = intersection(*Z_out, *comp);
            }
        }
        return Z_out;
    }

    std::unique_ptr<HybZono> ConZono::do_complement(const zono_float delta_m, bool, const OptSettings&, OptSolution*,
                                                    int, int)
    {
        // make sure in [-1,1] form
        if (this->is_0_1_form()) this->convert_form();

        // get a value lambda_m such that lambda_m >= max{ ||lambda||_infty : |[G^T A^T] lambda| <= 1 }

        // construct the matrix [G^T A^T]
        const auto GT = this->G.transpose();
        const auto AT = this->A.transpose();
        const Eigen::SparseMatrix<zono_float, Eigen::RowMajor> GTAT = hcat<zono_float>(GT, AT); // convert to row-major

        // get the smallest non-zero element in the matrix, making sure that all rows have at least one non-zero element
        zono_float min_non_zero = std::numeric_limits<zono_float>::max();
        for (int row = 0; row < GTAT.outerSize(); ++row)
        {
            bool value_exists = false;
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(GTAT, row); it; ++it)
            {
                if (it.value() < min_non_zero && std::abs(it.value()) > zono_eps)
                {
                    min_non_zero = it.value();
                }
                value_exists = true;
            }
            if (!value_exists)
            {
                std::stringstream ss;
                ss << "ConZono complement: row " << row << " of [G^T A^T] has no non-zero elements.";
                throw std::runtime_error(ss.str());
            }
        }

        // value for lambda_m
        const zono_float lambda_m = 1 / min_non_zero;

        // m value
        const zono_float m = delta_m + 1;

        // interval sets
        Eigen::Vector<zono_float, -1> l1(2 * this->nG);
        l1.setConstant(-(m + delta_m / 2));
        Eigen::Vector<zono_float, -1> u1(2 * this->nG);
        u1.setConstant(1 + delta_m / 2);
        auto Zf1 = interval_2_zono(Box(l1, u1));
        if (Zf1->is_0_1_form()) Zf1->convert_form();

        Eigen::Vector<zono_float, -1> l2(4 * this->nG);
        l2.segment(0, 2 * this->nG).setConstant(-(m + 3 * delta_m / 2 + 1));
        l2.segment(2 * this->nG, 2 * this->nG).setConstant(-2);
        Eigen::Vector<zono_float, -1> u2(4 * this->nG);
        u2.segment(0, 2 * this->nG).setConstant(delta_m / 2);
        u2.segment(2 * this->nG, 2 * this->nG).setZero();
        auto Zf2 = interval_2_zono(Box(l2, u2));
        if (Zf2->is_0_1_form()) Zf2->convert_form();

        // build complement

        // Gc
        Eigen::SparseMatrix<zono_float> Gc = m * this->G;
        Gc.conservativeResize(this->n, 9 * this->nG + this->n + this->nC + 1);

        // Gb
        Eigen::SparseMatrix<zono_float> Gb(this->n, 2 * this->nG);

        // c
        Eigen::Vector<zono_float, -1> c = this->c;

        // helper matrices

        // declarations
        std::vector<Eigen::Triplet<zono_float>> triplets;
        int n_offset = 0;

        // AcPF = [mI, -delta_m/2, 0, 0, 0;
        //         -mI, -delta_m/2, 0, 0, 0]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, this->nG, -delta_m / two);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, i, -m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG, -delta_m / two);
        }
        Eigen::SparseMatrix<zono_float> AcPF(2 * this->nG, 3 * this->nG + 1 + this->nC + this->n);
        AcPF.setFromTriplets(triplets.begin(), triplets.end());

        // AcDF = [0, 0, lambda_m [G^T A^T], 0.5 I, -0.5 I;
        //         0, 0, 0, 0.5 1^T, 0.5 1^T]
        triplets.clear();
        get_triplets_offset<zono_float>(lambda_m * GTAT, triplets, 0, this->nG + 1);
        n_offset = this->nG + 1 + this->n + this->nC;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, n_offset + i, p5);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG, n_offset + i, p5);
        }
        n_offset += this->nG;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, n_offset + i, -p5);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG, n_offset + i, p5);
        }
        Eigen::SparseMatrix<zono_float> AcDF(this->nG + 1, 3 * this->nG + 1 + this->nC + this->n);
        AcDF.setFromTriplets(triplets.begin(), triplets.end());

        // AcCS = [-mI, delta_m/2, 0, 0, 0;
        //         mI, delta_m/2, 0, 0, 0;
        //         0, 0, 0, I, 0;
        //         0, 0, 0, 0, I]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, -m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, this->nG, delta_m / two);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG, delta_m / two);
        }
        n_offset = this->nG + 1 + this->n + this->nC;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(2 * this->nG + i, n_offset + i, one);
        }
        n_offset += this->nG;
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(3 * this->nG + i, n_offset + i, one);
        }
        Eigen::SparseMatrix<zono_float> AcCS(4 * this->nG, 3 * this->nG + 1 + this->nC + this->n);
        AcCS.setFromTriplets(triplets.begin(), triplets.end());

        // AbCS = [mI, 0;
        //         0, mI;
        //         -I, 0;
        //         0, -I]
        triplets.clear();
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(i, i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(this->nG + i, this->nG + i, m);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(2 * this->nG + i, i, -one);
        }
        for (int i = 0; i < this->nG; ++i)
        {
            triplets.emplace_back(3 * this->nG + i, this->nG + i, -one);
        }
        Eigen::SparseMatrix<zono_float> AbCS(4 * this->nG, 2 * this->nG);
        AbCS.setFromTriplets(triplets.begin(), triplets.end());

        // bDF = [0;
        //        1-nG]
        Eigen::Vector<zono_float, -1> bDF(this->nG + 1);
        bDF.setZero();
        bDF(this->nG) = static_cast<zono_float>(1 - this->nG);

        // Ac = [mA, 0, 0, 0;
        //       AcPF, Gf1, 0;
        //       AcDF, 0, 0;
        //       AcCS, 0, Gf2]
        triplets.clear();
        int m_offset = 0;
        n_offset = this->nG + 1 + this->n + this->nC + 2 * this->nG;
        get_triplets_offset<zono_float>(m * this->A, triplets, 0, 0);
        m_offset += this->nC;
        get_triplets_offset<zono_float>(AcPF, triplets, m_offset, 0);
        get_triplets_offset<zono_float>(Zf1->G, triplets, m_offset, n_offset);
        m_offset += 2 * this->nG;
        get_triplets_offset<zono_float>(AcDF, triplets, m_offset, 0);
        m_offset += this->nG + 1;
        get_triplets_offset<zono_float>(AcCS, triplets, m_offset, 0);
        get_triplets_offset<zono_float>(Zf2->G, triplets, m_offset, n_offset + 2 * this->nG);
        Eigen::SparseMatrix<zono_float> Ac(7 * this->nG + this->nC + 1, 9 * this->nG + this->n + this->nC + 1);
        Ac.setFromTriplets(triplets.begin(), triplets.end());

        // Ab = [0;
        //       0;
        //       0;
        //       AbCS];
        triplets.clear();
        get_triplets_offset<zono_float>(AbCS, triplets, this->nC + 2 * this->nG + this->nG + 1, 0);
        Eigen::SparseMatrix<zono_float> Ab(7 * this->nG + this->nC + 1, 2 * this->nG);
        Ab.setFromTriplets(triplets.begin(), triplets.end());

        // b = [b;
        //      cf1;
        //      bDF;
        //      cf2]
        Eigen::Vector<zono_float, -1> b(7 * this->nG + this->nC + 1);
        b.segment(0, this->nC) = this->b;
        b.segment(this->nC, 2 * this->nG) = Zf1->c;
        b.segment(this->nC + 2 * this->nG, this->nG + 1) = bDF;
        b.segment(this->nC + 2 * this->nG + this->nG + 1, 4 * this->nG) = Zf2->c;

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, false, false);
    }

    std::unique_ptr<HybZono> set_diff(const HybZono& Z1, HybZono& Z2, const zono_float delta_m,
                                      const bool remove_redundancy,
                                      const OptSettings& settings, OptSolution* solution, const int n_leaves,
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

    // setup functions
    std::unique_ptr<HybZono> zono_union_2_hybzono(std::vector<Zono>& Zs, const bool expose_indicators)
    {
        // can't be empty
        if (Zs.empty())
        {
            throw std::invalid_argument("Zono union: empty input vector.");
        }

        // zonotope dimension
        int n_dims = Zs[0].n;
        const int n_zonos = static_cast<int>(Zs.size());

        // loop through Zs
        for (auto& Z : Zs)
        {
            // make sure dimensions are consistent
            if (Z.n != n_dims)
            {
                throw std::invalid_argument("Zono union: inconsistent dimensions.");
            }

            // convert to [0,1] form
            if (!Z.zero_one_form)
            {
                Z.convert_form();
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
            n_gens = Zs[i].nG;
            Gd = Zs[i].G.toDense();
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
                tripvec.emplace_back(j, i, Zs[i].c(j));
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

    std::unique_ptr<HybZono> vrep_2_hybzono(const std::vector<Eigen::Matrix<zono_float, -1, -1>>& Vpolys,
                                            const bool expose_indicators)
    {
        // error handling
        if (Vpolys.empty())
        {
            throw std::invalid_argument("set_from_vrep: Vpolys must have at least one polytope.");
        }

        // dimensions
        const int n_polys = static_cast<int>(Vpolys.size());
        const int n_dims = static_cast<int>(Vpolys[0].cols());
        int n_verts; // declare

        // check if all polytopes have the same number of dimensions
        for (const auto& Vpoly : Vpolys)
        {
            if (Vpoly.cols() != n_dims)
            {
                throw std::invalid_argument("set_from_vrep: all polytopes must have the same number of dimensions.");
            }
        }

        // initialize V and M matrices as std::vectors
        // each entry is a row
        std::vector<Eigen::Matrix<zono_float, 1, -1>> V_vec, M_vec;
        Eigen::Matrix<zono_float, 1, -1> M_row(n_polys);

        // loop through each polytope
        for (int i = 0; i < n_polys; i++)
        {
            n_verts = static_cast<int>(Vpolys[i].rows());
            for (int j = 0; j < n_verts; j++)
            {
                // check if the vertex is already in V_vec
                auto vertex_equal = [&](const Eigen::Matrix<zono_float, 1, -1>& v) -> bool
                {
                    return (v - Vpolys[i].row(j)).norm() < zono_eps;
                };

                if (auto it_V = std::find_if(V_vec.begin(), V_vec.end(), vertex_equal); it_V == V_vec.end())
                {
                    V_vec.emplace_back(Vpolys[i].row(j));
                    M_row.setZero();
                    M_row(i) = 1;
                    M_vec.push_back(M_row);
                }
                else
                {
                    const int idx = static_cast<int>(std::distance(V_vec.begin(), it_V));
                    M_vec[idx](i) = 1;
                }
            }
        }

        const int nV = static_cast<int>(V_vec.size()); // number of unique vertices

        // convert to Eigen matrices
        Eigen::Matrix<zono_float, -1, -1> V(n_dims, nV);
        Eigen::Matrix<zono_float, -1, -1> M(nV, n_polys);
        for (int i = 0; i < nV; i++)
        {
            V.col(i) = V_vec[i];
            M.row(i) = M_vec[i];
        }

        // directly build hybzono in [0,1] form

        // declare
        std::vector<Eigen::Triplet<zono_float>> tripvec;

        // output dimension
        int n_out = n_dims;
        if (expose_indicators)
            n_out += static_cast<int>(n_polys);

        // Gc = [V, 0]
        Eigen::SparseMatrix<zono_float> Gc = V.sparseView();
        Gc.conservativeResize(n_out, 2 * nV);

        // Gb = [0]
        tripvec.clear();
        for (int i = n_dims; i < n_out; ++i)
        {
            tripvec.emplace_back(i, i - n_dims, one);
        }
        Eigen::SparseMatrix<zono_float> Gb(n_out, n_polys);
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        Gb.setFromSortedTriplets(tripvec.begin(), tripvec.end());
#else
        Gb.setFromTriplets(tripvec.begin(), tripvec.end());
#endif

        // c = 0
        Eigen::Vector<zono_float, -1> c(n_out);
        c.setZero();

        // Ac = [1^T, 0^T;
        //       0^T, 0^T;
        //       I, diag[sum(M, 2)]]
        tripvec.clear();
        Eigen::SparseMatrix<zono_float> Ac(2 + nV, 2 * nV);
        Eigen::SparseMatrix<zono_float> I_nv(nV, nV);
        I_nv.setIdentity();
        for (int i = 0; i < nV; i++)
        {
            tripvec.emplace_back(0, i, one);
        }
        get_triplets_offset<zono_float>(I_nv, tripvec, 2, 0);
        Eigen::Vector<zono_float, -1> sum_M = M.rowwise().sum();
        for (int i = 0; i < nV; i++)
        {
            tripvec.emplace_back(2 + i, nV + i, sum_M(i));
        }
        Ac.setFromTriplets(tripvec.begin(), tripvec.end());

        // Ab = [0^T;
        //       1^T;
        //       -M]
        Eigen::SparseMatrix<zono_float> Ab(2 + nV, n_polys);
        tripvec.clear();
        for (int i = 0; i < n_polys; i++)
        {
            tripvec.emplace_back(1, i, one);
        }
        Eigen::SparseMatrix<zono_float> mM_sp = -M.sparseView();
        get_triplets_offset<zono_float>(mM_sp, tripvec, 2, 0);
        Ab.setFromTriplets(tripvec.begin(), tripvec.end());

        // b = [1;
        //      1;
        //      0]
        Eigen::Vector<zono_float, -1> b(2 + nV);
        b.setZero();
        b(0) = 1;
        b(1) = 1;

        // return hybrid zonotope
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, true, true);
    }

    std::unique_ptr<ConZono> vrep_2_conzono(const Eigen::Matrix<zono_float, -1, -1>& Vpoly)
    {
        // dimensions
        const int n_dims = static_cast<int>(Vpoly.cols());
        const int n_verts = static_cast<int>(Vpoly.rows());

        // make generators
        const Eigen::SparseMatrix<zono_float> G = Vpoly.transpose().sparseView();
        Eigen::Vector<zono_float, -1> c(n_dims);
        c.setZero();

        // make constraints
        std::vector<Eigen::Triplet<zono_float>> tripvec;
        Eigen::SparseMatrix<zono_float> A(1, n_verts);
        for (int i = 0; i < n_verts; i++)
        {
            tripvec.emplace_back(0, i, one);
        }
        A.setFromTriplets(tripvec.begin(), tripvec.end());

        Eigen::Vector<zono_float, -1> b(1);
        b(0) = 1;

        // return conzono
        return std::make_unique<ConZono>(G, c, A, b, true);
    }


    std::unique_ptr<Zono> interval_2_zono(const Box& box)
    {
        // generator matrix
        std::vector<Eigen::Triplet<zono_float>> triplets;
        Eigen::SparseMatrix<zono_float> G(static_cast<Eigen::Index>(box.size()), static_cast<Eigen::Index>(box.size()));
        for (int i = 0; i < static_cast<int>(box.size()); i++)
        {
            triplets.emplace_back(i, i, box[i].width() / two);
        }
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        G.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        G.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // center
        Eigen::Vector<zono_float, -1> c = box.center();

        // return zonotope
        return std::make_unique<Zono>(G, c, false);
    }


    std::unique_ptr<Zono> make_regular_zono_2D(const zono_float radius, const int n_sides, const bool outer_approx,
                                               const Eigen::Vector<zono_float, 2>& c)
    {
        // check number of sides
        if (n_sides % 2 != 0 || n_sides < 4)
        {
            throw std::invalid_argument("make_regular_zono_2D: number of sides must be even and >= 4.");
        }

        // check radius
        if (radius <= 0)
        {
            throw std::invalid_argument("make_regular_zono_2D: radius must be positive.");
        }

        // problem parameters
        const int n_gens = n_sides / 2;
        const zono_float dphi = pi / static_cast<zono_float>(n_gens);
        const zono_float R = outer_approx ? radius / std::cos(dphi / 2) : radius;

        // generator matrix
        const int n_gens_2 = n_gens / 2;
        zono_float phi = (static_cast<zono_float>(n_gens_2)) * dphi;
        const zono_float l_side = 2 * R * std::sin(dphi / 2);
        Eigen::Matrix<zono_float, -1, -1> G(2, n_gens);
        for (int i = 0; i < n_gens; i++)
        {
            G(0, i) = l_side * std::cos(phi);
            G(1, i) = l_side * std::sin(phi);
            phi -= dphi;
        }

        // return zonotope
        return std::make_unique<Zono>(p5 * G.sparseView(), c, false);
    }

    // convex relaxation
    std::unique_ptr<ConZono> HybZono::convex_relaxation() const
    {
        return std::make_unique<ConZono>(this->G, this->c, this->A, this->b, this->zero_one_form);
    }

    // bounding box
    Box HybZono::do_bounding_box(const OptSettings& settings, OptSolution* solution)
    {
        // if sharp, compute from convex relaxation
        if (this->sharp)
        {
            const auto Z_CR = this->convex_relaxation();
            return Z_CR->bounding_box(settings, solution);
        }

        // make sure dimension is at least 1
        if (this->n == 0)
        {
            throw std::invalid_argument("Bounding box: empty set");
        }

        // init search direction for bounding box
        Eigen::Vector<zono_float, -1> d(this->n);
        d.setZero();

        // declarations
        Box box(this->n); // init
        zono_float s_neg, s_pos;

        // build QP for ADMM
        const Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // mixed-integer solution

        // get support in all box directions
        for (int i = 0; i < this->n; i++)
        {
            // negative direction

            // update QP
            d.setZero();
            d(i) = -1;
            q = -this->G.transpose() * d;

            // solve
            OptSolution sol = this->mi_opt(P, q, 0, this->A, this->b, settings);
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_neg = -d.dot(this->G * sol.z + this->c);

            // positive direction

            // update QP
            d.setZero();
            d(i) = 1;
            q = -this->G.transpose() * d;

            // solve
            sol = this->mi_opt(P, q, 0, this->A, this->b, settings);
            if (sol.infeasible)
                throw std::invalid_argument("Bounding box: Z is empty");
            else
                s_pos = d.dot(this->G * sol.z + this->c);

            // store bounds
            box[i] = Interval(s_neg, s_pos);
        }

        return box;
    }

    zono_float HybZono::do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings& settings,
                                   OptSolution* solution)
    {
        // check dimensions
        if (this->n != d.size())
        {
            throw std::invalid_argument("Support: inconsistent dimensions.");
        }

        // if sharp, can solve as convex optimization problem
        if (this->sharp)
        {
            const auto Zc = this->convex_relaxation();
            return Zc->support(d, settings, solution);
        }

        // cost
        Eigen::SparseMatrix<zono_float> P(this->nG, this->nG);
        Eigen::Vector<zono_float, -1> q = -this->G.transpose() * d;

        // solve MIQP
        const OptSolution sol = this->mi_opt(P, q, 0, this->A, this->b, settings, solution);

        // check feasibility and return solution
        if (sol.infeasible) // Z is empty
            throw std::invalid_argument("Support: infeasible");
        else
            return d.dot(this->G * sol.z + this->c);
    }

    std::vector<ConZono> HybZono::get_leaves(const bool remove_redundancy, const OptSettings& settings,
                                             OptSolution* solution, const int n_leaves, const int contractor_iter) const
    {
        // allocate all threads to branch and bound
        OptSettings settings_get_leaves = settings;

        // get leaves as conzonos
        const std::vector<Eigen::Vector<zono_float, -1>> bin_leaves = this->get_bin_leaves(
            settings_get_leaves, solution, n_leaves);
        std::vector<ConZono> leaves;
        for (auto& xi_b : bin_leaves)
        {
            Eigen::Vector<zono_float, -1> cp = this->c + this->Gb * xi_b;
            Eigen::Vector<zono_float, -1> bp = this->b - this->Ab * xi_b;
            leaves.emplace_back(this->Gc, cp, this->Ac, bp, this->zero_one_form);
        }
        if (remove_redundancy)
        {
            for (auto& leaf : leaves)
            {
                leaf.remove_redundancy(contractor_iter);
            }
        }

        return leaves;
    }


    void ConZono::constraint_reduction()
    {
        // make sure there are constraints to remove
        if (this->nC == 0) return;

        // put set into [-1, 1] form
        if (this->zero_one_form) this->convert_form();

        // execute algorithm 1 from paper
        Eigen::Vector<zono_float, -1> x_lb(this->nG);
        Eigen::Vector<zono_float, -1> x_ub(this->nG);
        x_lb.setConstant(-1);
        x_ub.setConstant(1);
        Box E(x_lb, x_ub);
        x_lb.setConstant(-std::numeric_limits<zono_float>::infinity());
        x_ub.setConstant(std::numeric_limits<zono_float>::infinity());
        Box R(x_lb, x_ub);
        Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_rm = this->A;
        for (int i = 0; i < this->nC; ++i)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_j(A_rm, i); it_j; ++it_j)
            {
                const zono_float a_ij = it_j.value();
                Interval y(this->b(i) / a_ij, this->b(i) / a_ij);
                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_k(A_rm, i); it_k; ++it_k)
                {
                    if (it_j.col() == it_k.col()) continue;
                    y = y - E[it_k.col()].to_interval() * (it_k.value() / a_ij);
                }
                R[it_j.col()].intersect_assign(R[it_j.col()], y.as_view());
                E[it_j.col()].intersect_assign(E[it_j.col()], R[it_j.col()]);
            }
        }

        // make sure conzono isn't empty (interval check)
        for (int j = 0; j < this->nG; ++j)
        {
            if (E[j].is_empty())
                throw std::runtime_error("ConZono constraint reduction: set is empty");
        }

        // build Q matrix
        Eigen::SparseMatrix<zono_float> Q(this->nG + this->nC, this->nG + this->nC);

        Eigen::SparseMatrix<zono_float> I_nG(this->nG, this->nG);
        I_nG.setIdentity();
        const Eigen::SparseMatrix<zono_float> Phi = this->G.transpose() * this->G + I_nG;

        std::vector<Eigen::Triplet<zono_float>> triplets;
        get_triplets_offset<zono_float>(Phi, triplets, 0, 0);
        get_triplets_offset<zono_float>(this->A, triplets, this->nG, 0);
        get_triplets_offset<zono_float>(this->A.transpose(), triplets, 0, this->nG);
        Q.setFromTriplets(triplets.begin(), triplets.end());

        // factorize Q
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<zono_float>> Q_ldlt(Q);
        if (Q_ldlt.info() != Eigen::Success)
            throw std::runtime_error(
                "ConZono constraint reduction: Q matrix factorization failed, most likely A is not full row rank.");

        // get estimated Hausdorff error for eliminating each generator
        std::vector<std::pair<int, zono_float>> haus_vec; // (index, error)
        haus_vec.reserve(this->nG);
        auto shift_permute = [size=this->nG + this->nC + 1](const int start_index,
                                                            const int end_index) -> Eigen::PermutationMatrix<
            Eigen::Dynamic, Eigen::Dynamic>
        {
            assert(
                start_index >= 0 && end_index >= 0 && start_index < size && end_index < size && start_index <
                end_index);
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(size);
            for (int i = 0; i < start_index; ++i)
            {
                P.indices()[i] = i;
            }
            P.indices()[start_index] = end_index;
            for (int i = start_index + 1; i <= end_index; ++i)
            {
                P.indices()[i] = i - 1;
            }
            for (int i = end_index + 1; i < size; ++i)
            {
                P.indices()[i] = i;
            }
            return P;
        };
        Eigen::Vector<zono_float, -1> e_j(this->nG + this->nC);
        for (int j = 0; j < this->nG; ++j)
        {
            // get r_j
            const zono_float r_j = std::max<zono_float>(
                zero, std::max<zono_float>(std::abs(R[j].to_interval().lb), std::abs(R[j].to_interval().ub)) - one);
            if (r_j < zono_eps)
            {
                haus_vec.emplace_back(j, zero);
                continue;
            }

            // linear system matrix
            e_j.setZero();
            e_j(j) = 1;
            const Eigen::Vector<zono_float, -1> Qinv_e_j = Q_ldlt.solve(e_j);

            triplets.clear();
            for (int i = 0; i < this->nG + this->nC; ++i)
            {
                triplets.emplace_back(i, i, one);
                if (i == j)
                {
                    triplets.emplace_back(this->nG + this->nC, j, one); // extra row for e_j^T
                }
            }
            for (int i = 0; i < this->nG + this->nC; ++i)
            {
                triplets.emplace_back(i, this->nG + this->nC, Qinv_e_j(i));
            }
            Eigen::SparseMatrix<zono_float> M(this->nG + this->nC + 1, this->nG + this->nC + 1);
#if EIGEN_VERSION_AT_LEAST(5,0,0)
            M.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
            M.setFromTriplets(triplets.begin(), triplets.end());
#endif

            // RHS for linear system
            Eigen::Vector<zono_float, -1> rhs(this->nG + this->nC + 1);
            rhs.setZero();
            rhs(this->nG + this->nC) = r_j;

            // permutation matrices to make system upper triangular
            const auto P_R = shift_permute(j, this->nG + this->nC);
            const auto P_L = shift_permute(j, this->nG + this->nC - 1);

            // solve linear system in permuted space
            const Eigen::SparseMatrix<zono_float> M_perm = P_L * M * P_R.inverse();
            const Eigen::Vector<zono_float, -1> rhs_perm = P_L * rhs;
            const Eigen::Vector<zono_float, -1> y_perm = M_perm.triangularView<Eigen::Upper>().solve(rhs_perm);
            const Eigen::Vector<zono_float, -1> y = P_R * y_perm;

            // get Hausdorff distance estimate
            const Eigen::Vector<zono_float, -1> d = y.segment(0, this->nG);
            const zono_float haus_j = (this->G * d).norm() + d.norm();
            haus_vec.emplace_back(j, haus_j);
        }

        // sort by Hausdorff error
        std::sort(haus_vec.begin(), haus_vec.end(), [](const auto& a, const auto& b) { return a.second < b.second; });

        Eigen::SparseMatrix<zono_float> Ea(this->nG, this->nC); // init
        int gen_remove = -1;
        int cons_remove = -1;
        for (const auto& [j, err] : haus_vec)
        {
            // loop through column k of A to find a constraint to remove
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->A, j); it; ++it)
            {
                if (std::abs(it.value()) > zono_eps)
                {
                    triplets.emplace_back(j, static_cast<int>(it.row()), one / it.value());
                    Ea.insert(j, it.row()) = 1 / it.value();
                    gen_remove = j;
                    cons_remove = static_cast<int>(it.row());
                    break;
                }
            }

            // check if done
            if (gen_remove != -1 && cons_remove != -1)
                break;
        }
        if (gen_remove == -1 || cons_remove == -1)
            throw std::runtime_error("ConZono: constraint reduction cannot find valid constraint to remove");

        // apply algorithm from Scott paper
        const Eigen::SparseMatrix<zono_float> Lambda_G = G * Ea;
        const Eigen::SparseMatrix<zono_float> Lambda_A = A * Ea;
        Eigen::SparseMatrix<zono_float> Gp = this->G - Lambda_G * this->A;
        Eigen::Vector<zono_float, -1> cp = this->c + Lambda_G * this->b;
        Eigen::SparseMatrix<zono_float> Ap = this->A - Lambda_A * this->A;
        Eigen::Vector<zono_float, -1> bp = this->b - Lambda_A * this->b;

        // generator removal matrix
        triplets.clear();
        for (int j = 0; j < gen_remove; ++j)
        {
            triplets.emplace_back(j, j, one);
        }
        for (int j = gen_remove + 1; j < this->nG; ++j)
        {
            triplets.emplace_back(j, j - 1, one);
        }
        Eigen::SparseMatrix<zono_float> dG(this->nG, this->nG - 1);
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        dG.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        dG.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // constraint removal matrix
        triplets.clear();
        for (int i = 0; i < cons_remove; ++i)
        {
            triplets.emplace_back(i, i, one);
        }
        for (int i = cons_remove + 1; i < this->nC; ++i)
        {
            triplets.emplace_back(i - 1, i, one);
        }
        Eigen::SparseMatrix<zono_float> dA(this->nC - 1, this->nC);
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        dA.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        dA.setFromTriplets(triplets.begin(), triplets.end());
#endif

        // update
        this->set(Gp * dG, cp, dA * Ap * dG, dA * bp, false);
    }

    std::unique_ptr<Zono> ConZono::to_zono_approx() const
    {
        // check for case that there are no constraints
        if (this->nG == 0)
        {
            return std::make_unique<Point>(this->c);
        }
        if (this->nC == 0)
        {
            return std::make_unique<Zono>(this->G, this->c, this->zero_one_form);
        }

        // compute SVD of A
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        const Eigen::BDCSVD<Eigen::Matrix<zono_float, -1, -1>, Eigen::ComputeFullV | Eigen::ComputeFullU> svd(
            this->A.toDense());
#else
        const Eigen::BDCSVD<Eigen::Matrix<zono_float, -1, -1>> svd(
            this->A.toDense(), Eigen::ComputeFullV | Eigen::ComputeFullU);
#endif
        const Eigen::Vector<zono_float, -1>& sin_vals = svd.singularValues();

        const int n = this->nG;
        const int r = static_cast<int>(svd.rank());
        const int d = n - r;
        const Eigen::Matrix<zono_float, -1, -1> Vr = svd.matrixV().block(0, 0, n, r);
        const Eigen::Matrix<zono_float, -1, -1> Vd = svd.matrixV().block(0, r, n, d);

        // get xi_tilde_r
        const Eigen::Vector<zono_float, -1> b_tilde = svd.matrixU().transpose() * this->b;
        const Eigen::Vector<zono_float, -1> xi_tilde_r = b_tilde.segment(0, r).array() / sin_vals.segment(0, r).array();

        // build bounding zonotope
        const Eigen::SparseMatrix<zono_float> Vd_VdT_sp = (Vd * Vd.transpose()).sparseView();
        const Eigen::SparseMatrix<zono_float> G_zono = this->G * Vd_VdT_sp;
        const Eigen::Vector<zono_float, -1> c_zono = this->c + this->G * (Vr * xi_tilde_r);
        return std::make_unique<Zono>(G_zono, c_zono, this->zero_one_form);
    }
} // namespace ZonoOpt
