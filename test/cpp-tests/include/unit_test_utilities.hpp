#ifndef UNIT_TEST_UTILITIES_HPP_
#define UNIT_TEST_UTILITIES_HPP_

#include "ZonoOpt.hpp"
#include <fstream>
#include <iostream>
#include <random>

inline Eigen::SparseMatrix<zono_float> load_sparse_matrix(const std::string& filename)
{
    std::ifstream file(filename);
    std::string line;
    int row = 0, col = 0;
    std::vector<Eigen::Triplet<zono_float>> triplets;

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        double val;
        col = 0;
        while (ss >> val)
        {
            if (std::abs(val) > zono_eps)
                triplets.emplace_back(row, col, static_cast<zono_float>(val));
            ++col;
        }
        ++row;
    }

    Eigen::SparseMatrix<zono_float> mat(row, col);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

inline Eigen::Vector<zono_float, -1> load_vector(const std::string& filename)
{
    std::ifstream file(filename);
    std::string line;
    std::vector<zono_float> values;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        double val;
        while (ss >> val)
        {
            values.push_back(static_cast<zono_float>(val));
        }
    }
    Eigen::Vector<zono_float, -1> vec(values.size());
    for (size_t i = 0; i < values.size(); ++i)
    {
        vec(i) = values[i];
    }
    return vec;
}

inline void test_assert(const bool assert_cond, const std::string& message = "")
{
    if (!assert_cond)
    {
        std::cerr << "Assertion failed: " <<  message << std::endl;
        std::exit(1);
    }
}

inline Eigen::SparseMatrix<zono_float> random_sparse_matrix(int m, int n, double density, zono_float val_min, zono_float val_max, std::mt19937& rand_gen)
{
    // input checking
    if (m < 0 || n < 0)
        throw std::invalid_argument("matrix dimensions must be non-negative");
    if (density <= 0 || density > 1)
        throw std::invalid_argument("density must be in (0, 1]");
    if (val_max < val_min)
        throw std::invalid_argument("val_max must be greater than val_min");

    // initialize randomization
    std::uniform_real_distribution<zono_float> val_dist(val_min, val_max);
    std::uniform_real_distribution<double> density_dist(0.0, 1.0);

    // generate random triplets
    Eigen::SparseMatrix<zono_float> mat(m, n);
    std::vector<Eigen::Triplet<zono_float>> triplets;
    for (int i=0; i<m; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            if (density_dist(rand_gen) < density)
            {
                triplets.emplace_back(i, j, val_dist(rand_gen));
            }
        }
    }
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
    mat.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
    mat.setFromTriplets(triplets.begin(), triplets.end());
#endif

    return mat;
}

inline Eigen::Vector<zono_float, -1> random_vector(int n, zono_float val_min, zono_float val_max, std::mt19937& rand_gen)
{
    // input checking
    if (n < 0)
        throw std::invalid_argument("vector dimension must be non-negative");
    if (val_max < val_min)
        throw std::invalid_argument("val_max must be greater than val_min");

    // initialize randomization
    std::uniform_real_distribution<zono_float> val_dist(val_min, val_max);

    // generate random vector
    Eigen::Vector<zono_float, -1> vec(n);
    for (int i=0; i<n; ++i)
    {
        vec(i) = val_dist(rand_gen);
    }

    return vec;
}

inline ZonoOpt::HybZono random_hybzono(int n, int nGc, int nGb, int nC, double density, zono_float val_min, zono_float val_max, std::mt19937& rand_gen)
{
    Eigen::SparseMatrix<zono_float> Gc = random_sparse_matrix(n, nGc, density, val_min, val_max, rand_gen);
    Eigen::SparseMatrix<zono_float> Gb = random_sparse_matrix(n, nGb, density, val_min, val_max, rand_gen);
    Eigen::Vector<zono_float, -1> c = random_vector(n, val_min, val_max, rand_gen);
    Eigen::SparseMatrix<zono_float> Ac = random_sparse_matrix(nC, nGc, density, val_min, val_max, rand_gen);
    Eigen::SparseMatrix<zono_float> Ab = random_sparse_matrix(nC, nGb, density, val_min, val_max, rand_gen);
    Eigen::Vector<zono_float, -1> b = random_vector(nC, val_min, val_max, rand_gen);

    return {Gc, Gb, c, Ac, Ab, b};
}

inline ZonoOpt::ConZono random_conzono(int n, int nG, int nC, double density, zono_float val_min, zono_float val_max, std::mt19937& rand_gen)
{
    Eigen::SparseMatrix<zono_float> G = random_sparse_matrix(n, nG, density, val_min, val_max, rand_gen);
    Eigen::Vector<zono_float, -1> c = random_vector(n, val_min, val_max, rand_gen);
    Eigen::SparseMatrix<zono_float> A = random_sparse_matrix(nC, nG, density, val_min, val_max, rand_gen);
    Eigen::Vector<zono_float, -1> b = random_vector(nC, val_min, val_max, rand_gen);

    return {G, c, A, b};
}

inline ZonoOpt::Zono random_zono(int n, int nG, double density, zono_float val_min, zono_float val_max, std::mt19937& rand_gen)
{
    Eigen::SparseMatrix<zono_float> G = random_sparse_matrix(n, nG, density, val_min, val_max, rand_gen);
    Eigen::Vector<zono_float, -1> c = random_vector(n, val_min, val_max, rand_gen);

    return {G, c};
}

inline bool matrix_eq(const Eigen::Matrix<zono_float, -1, -1>& A, const Eigen::Matrix<zono_float, -1, -1>& B)
{
    if (A.rows() != B.rows() || A.cols() != B.cols())
        return false;

    for (int i=0; i<A.rows(); ++i)
    {
        for (int j=0; j<A.cols(); ++j)
        {
            if (std::abs(A(i,j) - B(i,j)) > zono_eps)
                return false;
        }
    }
    return true;
}

inline bool matrix_eq(const Eigen::SparseMatrix<zono_float>& A, const Eigen::SparseMatrix<zono_float>& B)
{
    // convert to dense
    const Eigen::Matrix<zono_float, -1, -1>& Ad = A.toDense();
    const Eigen::Matrix<zono_float, -1, -1>& Bd = B.toDense();
    return matrix_eq(Ad, Bd);
}

inline bool vector_eq(const Eigen::Vector<zono_float, -1>& A, const Eigen::Vector<zono_float, -1>& B)
{
    if (A.size() != B.size())
        return false;

    for (int i=0; i<A.size(); ++i)
    {
        if (std::abs(A(i) - B(i)) > zono_eps)
            return false;
    }
    return true;
}

inline bool hz_eq(const ZonoOpt::HybZono& A, const ZonoOpt::HybZono& B)
{
    return matrix_eq(A.get_Gc(), B.get_Gc()) && matrix_eq(A.get_Gb(), B.get_Gb()) && vector_eq(A.get_c(), B.get_c())
        && matrix_eq(A.get_Ac(), B.get_Ac()) && matrix_eq(A.get_Ab(), B.get_Ab()) && vector_eq(A.get_b(), B.get_b())
        && A.is_0_1_form() == B.is_0_1_form();
}


#endif