#ifndef ZONOOPT_ZONO_JSON_HPP__
#define ZONOOPT_ZONO_JSON_HPP__

/**
 * @file ZonoJson.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Read / write zonotopic sets from / to JSON files
 * @version 1.0
 * @date 2026-04-12
 *
 * @copyright Copyright (c) 2026
 *
 */


#include <fstream>
#include <string>
#include "nlohmann/json.hpp"
#include "HybZono.hpp"

using json = nlohmann::json;

namespace ZonoOpt
{
    /**
     * @brief Write a HybZono object to a JSON file
     *
     * @param Z zonotopic set
     * @param filename the name of the JSON file to write to
     * @ingroup ZonoOpt_SetupFunctions
     */
    void to_json(const HybZono& Z, const std::string& filename);

    /**
     * @brief Read a HybZono object from a JSON file
     *
     * @param filename The name of the JSON file to read from
     * @return zonotopic set
     * @ingroup ZonoOpt_SetupFunctions
     */
    std::unique_ptr<HybZono> from_json(const std::string& filename);

namespace detail
{
    void sparse_to_json(const Eigen::SparseMatrix<zono_float>& mat, json& j);
    Eigen::SparseMatrix<zono_float> json_to_sparse(const json& j);
    Eigen::Vector<zono_float, -1> std_2_eigen(const std::vector<zono_float>& vec);
    std::vector<zono_float> eigen_2_std(const Eigen::Vector<zono_float, -1>& vec);
}

}


#endif