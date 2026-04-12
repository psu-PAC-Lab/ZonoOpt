#include "ZonoOpt.hpp"

namespace ZonoOpt
{
using namespace detail;

void to_json(const HybZono& Z, const std::string& filename)
{
    // open file to write to
    std::ofstream file(filename.c_str());

    // create json object
    json data;

    // class string
    std::string class_str;
    if (Z.is_hybzono())
        class_str = "HybZono";
    else if (Z.is_conzono())
        class_str = "ConZono";
    else if (Z.is_zono())
        class_str = "Zono";
    else if (Z.is_point())
        class_str = "Point";
    else if (Z.is_empty_set())
        class_str = "EmptySet";
    else
        throw std::runtime_error("to_json: unrecognized set type.");

    // serialize data
    data["class"] = class_str;
    data["n"] = Z.get_n();
    data["Gc"] = json();
    sparse_to_json(Z.get_Gc(), data["Gc"]);
    data["Gb"] = json();
    sparse_to_json(Z.get_Gb(), data["Gb"]);
    if (Z.is_empty_set())
        data["c"] = std::vector<zono_float>(Z.get_n(), zero); // handle NaNs in c for EmptySet
    else
        data["c"] = eigen_2_std(Z.get_c());
    data["Ac"] = json();
    sparse_to_json(Z.get_Ac(), data["Ac"]);
    data["Ab"] = json();
    sparse_to_json(Z.get_Ab(), data["Ab"]);
    data["b"] = eigen_2_std(Z.get_b());
    data["zero_one_form"] = Z.is_0_1_form();

    // write to file
    if (file)
        file << data;
    else
        throw std::runtime_error("to_json: could not open file for writing.");
}

std::unique_ptr<HybZono> from_json(const std::string& filename)
{
    // open file to read from
    std::ifstream file(filename.c_str());

    // read to json object
    json data;
    if (file)
        file >> data;
    else
        throw std::runtime_error("from_json: could not open file for reading.");

    // deserialize data
    const Eigen::SparseMatrix<zono_float> Gc = json_to_sparse(data["Gc"]);
    const Eigen::SparseMatrix<zono_float> Gb = json_to_sparse(data["Gb"]);
    const Eigen::Vector<zono_float, -1> c = std_2_eigen(data["c"]);
    const Eigen::SparseMatrix<zono_float> Ac = json_to_sparse(data["Ac"]);
    const Eigen::SparseMatrix<zono_float> Ab = json_to_sparse(data["Ab"]);
    const Eigen::Vector<zono_float, -1> b = std_2_eigen(data["b"]);
    const bool zero_one_form = data["zero_one_form"];
    const int n = data["n"];

    // return zonotopic set
    std::string class_str = data["class"];
    std::string hz_class_str = "HybZono";
    if (data["class"] == "HybZono")
        return std::make_unique<HybZono>(Gc, Gb, c, Ac, Ab, b, zero_one_form);
    else if (data["class"] == "ConZono")
        return std::make_unique<ConZono>(Gc, c, Ac, b, zero_one_form);
    else if (data["class"] == "Zono")
        return std::make_unique<Zono>(Gc, c, zero_one_form);
    else if (data["class"] == "Point")
        return std::make_unique<Point>(c);
    else if (data["class"] == "EmptySet")
        return std::make_unique<EmptySet>(n);
    else
        throw std::runtime_error("from_json: unrecognized set type.");
}

namespace detail
{
    void sparse_to_json(const Eigen::SparseMatrix<zono_float>& mat, json& j)
    {
        j["rows"] = mat.rows();
        j["cols"] = mat.cols();

        std::vector<int> trip_rows(mat.nonZeros());
        std::vector<int> trip_cols(mat.nonZeros());
        std::vector<zono_float> trip_vals(mat.nonZeros());

        int idx = 0;
        for (int k=0; k<mat.outerSize(); ++k)
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(mat,k); it; ++it)
            {
                trip_vals[idx] = it.value();
                trip_rows[idx] = static_cast<int>(it.row());
                trip_cols[idx] = static_cast<int>(it.col());
                ++idx;
            }

        j["trip_rows"] = trip_rows;
        j["trip_cols"] = trip_cols;
        j["trip_vals"] = trip_vals;
    }

    Eigen::SparseMatrix<zono_float> json_to_sparse(const json& j)
    {
        const int rows = j["rows"];
        const int cols = j["cols"];
        const std::vector<int> trip_rows(j["trip_rows"]);
        const std::vector<int> trip_cols(j["trip_cols"]);
        const std::vector<zono_float> trip_vals(j["trip_vals"]);

        if (trip_rows.size() != trip_cols.size() || trip_rows.size() != trip_vals.size())
            throw std::invalid_argument("json_to_sparse: triplet vectors must have the same size.");

        Eigen::SparseMatrix<zono_float> mat (rows, cols);
        std::vector<Eigen::Triplet<zono_float>> triplets;
        triplets.reserve(trip_rows.size());
        for (int i=0; i<static_cast<int>(trip_rows.size()); ++i)
        {
            triplets.emplace_back(trip_rows[i], trip_cols[i], trip_vals[i]);
        }
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        mat.setFromSortedTriplets(triplets.begin(), triplets.end());
#else
        mat.setFromTriplets(triplets.begin(), triplets.end());
#endif

        return mat;
    }

    Eigen::Vector<zono_float, -1> std_2_eigen(const std::vector<zono_float>& vec)
    {
        Eigen::Vector<zono_float, -1> v_eig (vec.size());
        for (int i=0; i<static_cast<int>(vec.size()); ++i)
        {
            v_eig[i] = vec[i];
        }
        return v_eig;
    }

    std::vector<zono_float> eigen_2_std(const Eigen::Vector<zono_float, -1>& vec)
    {
        std::vector<zono_float> v_std (vec.size());
        for (int i=0; i<static_cast<int>(vec.size()); ++i)
        {
            v_std[i] = vec[i];
        }
        return v_std;
    }

}

}