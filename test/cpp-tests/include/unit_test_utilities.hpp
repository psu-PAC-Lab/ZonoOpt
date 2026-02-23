#ifndef UNIT_TEST_UTILITIES_HPP_
#define UNIT_TEST_UTILITIES_HPP_

#include <fstream>
#include "ZonoOpt.hpp"

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



#endif