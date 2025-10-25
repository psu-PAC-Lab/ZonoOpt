#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    Zono::Zono(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
               const bool zero_one_form)
    {
        set(G, c, zero_one_form);
        sharp = true;
    }

    HybZono* Zono::clone() const
    {
        return new Zono(*this);
    }

    void Zono::set(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                   const bool zero_one_form)
    {
        // check dimensions
        if (G.rows() != c.size())
        {
            throw std::invalid_argument("Zono: inconsistent dimensions.");
        }

        // zonotope parameters
        this->G = G;
        this->c = c;
        this->nG = static_cast<int>(this->G.cols());
        this->n = static_cast<int>(this->G.rows());
        this->zero_one_form = zero_one_form;

        // abstract zono parameters
        this->nGc = this->nG;
        this->nGb = 0;
        this->nC = 0;
        this->Gc = this->G;
        this->Gb.resize(this->n, 0);
        this->A.resize(0, this->nG);
        this->Ac = this->A;
        this->Ab.resize(0, 0);
        this->b.resize(0);
    }

    void Zono::convert_form()
    {
        Eigen::Vector<zono_float, -1> c;
        Eigen::SparseMatrix<zono_float> G;

        if (!this->zero_one_form) // convert to [0,1] generators
        {
            c = this->c - this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            G = 2.0 * this->G;

            set(G, c, true);
        }
        else // convert to [-1,1] generators
        {
            c = this->c + 0.5 * this->G * Eigen::Vector<zono_float, -1>::Ones(this->nG);
            G = 0.5 * this->G;

            set(G, c, false);
        }
    }

    std::string Zono::print() const
    {
        std::stringstream ss;
        ss << "Zono: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "nG: " << this->nG << std::endl;
        ss << "G: " << Eigen::Matrix<zono_float, -1, -1>(this->G) << std::endl;
        ss << "c: " << this->c << std::endl;
        ss << "zero_one_form: " << this->zero_one_form;
        return ss.str();
    }

    bool Zono::do_is_empty(const OptSettings&, OptSolution*) const
    {
        if (this->n == 0)
            return true;
        else
            return false;
    }

    Box Zono::do_bounding_box(const OptSettings&, OptSolution*)
    {
        // convert to [-1,1] form
        if (this->zero_one_form) this->convert_form();

        // init bounds
        Eigen::Vector<zono_float, -1> l = this->c;
        Eigen::Vector<zono_float, -1> u = this->c;

        // compute bounds
        for (int i = 0; i < this->nG; ++i)
        {
            l -= this->G.col(i).cwiseAbs();
            u += this->G.col(i).cwiseAbs();
        }

        // return zonotope bounding box
        return {l, u};
    }

    std::unique_ptr<Zono> Zono::reduce_order(const int n_o)
    {
        // check validity
        if (n_o < this->n)
            throw std::invalid_argument("Zono reduce_order: desired order is less than dimension of set");

        // trivial case
        if (this->nG <= n_o)
            return std::make_unique<Zono>(*this);

        // convert to [-1,1] form
        if (this->zero_one_form) this->convert_form();

        // sort columns by decreasing 2-norm

        // vector of 2-norms and indices
        std::vector<std::pair<int, zono_float>> sort_vec;
        for (int i = 0; i < this->nG; ++i)
        {
            sort_vec.emplace_back(i, this->G.col(i).norm());
        }

        // sort
        auto comp = [](const std::pair<int, zono_float>& a, const std::pair<int, zono_float>& b) -> bool
        {
            return a.second > b.second;
        };
        std::sort(sort_vec.begin(), sort_vec.end(), comp);

        // zonotope to keep
        const int n_K = n_o - this->n;
        Eigen::SparseMatrix<zono_float> G_K(this->n, n_K);
        std::vector<Eigen::Triplet<zono_float>> triplets;
        for (int i = 0; i < n_K; ++i)
        {
            const int k = sort_vec[static_cast<size_t>(i)].first; // column
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->G, k); it; ++it)
            {
                triplets.emplace_back(static_cast<int>(it.row()), i, it.value());
            }
        }
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        G_K.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        G_K.setFromTriplets(triplets.begin(), triplets.end());
#endif
        const Zono K(G_K, this->c);

        // zonotope to over-approximate
        Eigen::SparseMatrix<zono_float> G_L(this->n, this->nG - n_K);
        triplets.clear();
        for (int i = n_K; i < this->nG; ++i)
        {
            const int k = sort_vec[static_cast<size_t>(i)].first; // column
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(this->G, k); it; ++it)
            {
                triplets.emplace_back(static_cast<int>(it.row()), i - n_K, it.value());
            }
        }
#if EIGEN_VERSION_AT_LEAST(5,0,0)
        G_L.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        G_L.setFromTriplets(triplets.begin(), triplets.end());
#endif
        Zono L(G_L, Eigen::Vector<zono_float, -1>::Zero(this->n));

        // get bounding box
        const auto L_R = interval_2_zono(L.bounding_box());

        // minkowski sum
        auto Z = minkowski_sum(K, *L_R);

        // check that dynamic cast is valid
        if (!Z->is_zono())
            throw std::runtime_error("In Zono::ReduceOrder, return type is not a zonotope?");

        // cast to zono and return
        return std::unique_ptr<Zono>(dynamic_cast<Zono*>(Z.release()));
    }

    zono_float Zono::do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings&, OptSolution*)
    {
        if (this->zero_one_form) this->convert_form();

        zono_float h = d.dot(this->c);
        const Eigen::Matrix<zono_float, -1, -1> Gd = this->G.toDense();
        for (int i = 0; i < this->nG; ++i)
        {
            h += std::abs(d.dot(Gd.col(i)));
        }
        return h;
    }


    zono_float Zono::get_volume()
    {
        // require [0,1] form
        if (!this->zero_one_form) this->convert_form();

        // get dense form of generator matrix
        const Eigen::Matrix<zono_float, -1, -1> Gd = this->G.toDense();

        // get cols of generator matrix
        std::vector<int> cols;
        cols.reserve(static_cast<size_t>(this->nG));
        for (int i = 0; i < this->nG; ++i)
        {
            cols.push_back(i);
        }

        // get combinations
        const auto combs = get_combinations<int>(cols, static_cast<size_t>(this->n));

        // get determinants
        zono_float det_sum = zero;
        Eigen::LDLT<Eigen::Matrix<zono_float, -1, -1>> ldlt;

        for (const auto& comb : combs)
        {
#if EIGEN_VERSION_AT_LEAST(5,0,0)
            const Eigen::Matrix<zono_float, -1, -1> G_comb = Gd(Eigen::placeholders::all, comb);
#else
            const Eigen::Matrix<zono_float, -1, -1> G_comb = Gd(Eigen::all, comb);
#endif
            const Eigen::Matrix<zono_float, -1, -1> GT_G = G_comb.transpose() * G_comb;
            ldlt.compute(GT_G);
            const Eigen::Vector<zono_float, -1> D = ldlt.vectorD();
            det_sum += std::sqrt(D.prod());
        }

        // volume is det_sum
        return det_sum;
    }
}
