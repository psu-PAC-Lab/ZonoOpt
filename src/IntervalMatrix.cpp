#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

     IntervalMatrix::IntervalMatrix(const size_t rows, const size_t cols, const std::vector<Eigen::Triplet<Interval>>& triplets)
    {
        // store dimensions
        this->rows_ = rows;
        this->cols_ = cols;

        // build internal storage
        this->mat_.resize(rows);
        for (const auto& triplet : triplets)
        {
            auto it = std::find_if(this->mat_[triplet.row()].begin(), this->mat_[triplet.row()].end(),
                                   [&](const std::pair<size_t, Interval>& p)
                                   {
                                       return p.first == static_cast<size_t>(triplet.col());
                                   });
            if (it != this->mat_[triplet.row()].end())
            {
                it->second = it->second + triplet.value();
            }
            else
            {
                this->mat_[triplet.row()].emplace_back(triplet.col(), triplet.value());
            }
        }

        // sort each row by column index
        auto comp = [](const std::pair<size_t, Interval>& a, const std::pair<size_t, Interval>& b) -> bool
        {
            return a.first < b.first;
        };
        for (auto& row : this->mat_)
        {
            std::sort(row.begin(), row.end(), comp);
        }
    }

    IntervalMatrix::IntervalMatrix(const Eigen::Matrix<zono_float, -1, -1>& mat_lb,
                                   const Eigen::Matrix<zono_float, -1, -1>& mat_ub)
    {
        if (mat_lb.rows() != mat_ub.rows() || mat_lb.cols() != mat_ub.cols())
            throw std::invalid_argument("IntervalMatrix: mat_lb and mat_ub must have the same dimensions");

        // initialize rows / cols
        this->rows_ = mat_lb.rows();
        this->cols_ = mat_lb.cols();

        // get triplets
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < mat_lb.rows(); ++i)
        {
            for (int j = 0; j < mat_lb.cols(); j++)
            {
                if (std::abs(mat_lb(i, j)) > zono_eps || std::abs(mat_ub(i, j)) > zono_eps)
                    triplets.emplace_back(i, j, Interval(mat_lb(i, j), mat_ub(i, j)));
            }
        }
        *this = IntervalMatrix(static_cast<size_t>(mat_lb.rows()), static_cast<size_t>(mat_lb.cols()), triplets);
    }

    IntervalMatrix::IntervalMatrix(const Eigen::Matrix<Interval, -1, -1>& mat)
    {
        // initialize rows / cols
        this->rows_ = mat.rows();
        this->cols_ = mat.cols();

        // get triplets
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < mat.rows(); ++i)
        {
            for (int j = 0; j < mat.cols(); j++)
            {
                if (std::abs(mat(i, j).lower()) > zono_eps || std::abs(mat(i, j).upper()) > zono_eps)
                    triplets.emplace_back(i, j, mat(i, j));
            }
        }
        *this = IntervalMatrix(static_cast<size_t>(mat.rows()), static_cast<size_t>(mat.cols()), triplets);
    }

    Box IntervalMatrix::operator*(const Eigen::Vector<zono_float, -1>& v) const
    {
        // input handling
        if (v.size() != static_cast<Eigen::Index>(this->cols_))
            throw std::invalid_argument("IntervalMatrix multiplication with box: inconsistent dimensions");

        // declare
        Box y(this->rows_);

        // linear map
        for (int i = 0; i < static_cast<int>(this->rows_); i++)
        {
            y.set_element(i, Interval(0, 0));
            for (const auto& [j, val] : this->mat_[i])
            {
                y.set_element(i, y.get_element(i) + val * v[static_cast<Eigen::Index>(j)]);
            }
        }
        return y;
    }

    Box IntervalMatrix::operator*(const Box& box) const
    {
        // input handling
        if (box.size() != this->cols_)
            throw std::invalid_argument("IntervalMatrix multiplication with box: inconsistent dimensions");

        // declare
        Box y(this->rows_);

        // linear map
        for (int i = 0; i < static_cast<int>(this->rows_); i++)
        {
            y.set_element(i, Interval(0, 0));
            for (const auto& [j, val] : this->mat_[i])
            {
                y.set_element(i, y.get_element(i) + (box.get_element(static_cast<int>(j)) * val));
            }
        }
        return y;
    }

    IntervalMatrix operator+(const Interval& interval, const IntervalMatrix& mat)
    {
        return mat + interval;
    }

    IntervalMatrix IntervalMatrix::operator*(zono_float alpha) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val *= alpha;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator*=(zono_float alpha)
    {
        *this = *this * alpha;
    }

    IntervalMatrix IntervalMatrix::operator*(const Interval& interval) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val *= interval;
            }
        }
        return mat;
    }

    IntervalMatrix operator*(const Interval& interval, const IntervalMatrix& A)
    {
        return A * interval;
    }

    void IntervalMatrix::operator*=(const Interval& interval)
    {
        *this = *this * interval;
    }

    IntervalMatrix operator*(zono_float alpha, const IntervalMatrix& A)
    {
        return A * alpha;
    }

    IntervalMatrix IntervalMatrix::operator/(zono_float alpha) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val /= alpha;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator/=(zono_float alpha)
    {
        *this = *this / alpha;
    }

    IntervalMatrix IntervalMatrix::operator/(const Interval& interval) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val /= interval;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator/=(const Interval& interval)
    {
        *this = *this / interval;
    }

    IntervalMatrix IntervalMatrix::operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const
    {
        // input handling
        if (A.rows() != static_cast<Eigen::Index>(this->cols_))
            throw std::invalid_argument("IntervalMatrix multiplication with matrix: inconsistent dimensions");

        // dimensions
        size_t rows = this->rows_;
        auto cols = static_cast<size_t>(A.cols());

        // loop through to generate triplets
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < static_cast<int>(rows); ++i)
        {
            for (const auto& [k, val] : this->mat_[i])
            {
                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_A(
                         A, static_cast<Eigen::Index>(k)); it_A; ++it_A)
                {
                    triplets.emplace_back(i, static_cast<int>(it_A.col()), val * it_A.value());
                }
            }
        }
        return {rows, cols, triplets};
    }

    IntervalMatrix IntervalMatrix::operator*(const Eigen::Matrix<zono_float, -1, -1>& A) const
    {
        const Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_sp = A.sparseView();
        return *this * A_sp;
    }

    IntervalMatrix operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const IntervalMatrix& B) {
         // input handling
         if (A.cols() != static_cast<Eigen::Index>(B.rows()))
             throw std::invalid_argument("IntervalMatrix multiplication with matrix: inconsistent dimensions");

         // dimensions
         size_t rows = A.rows();
         size_t cols = B.cols();

         // loop through to generate triplets
         std::vector<Eigen::Triplet<Interval>> triplets;
         for (int i=0; i<A.rows(); ++i) {
             for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it) {
                for (const auto& [j, val] : B.mat_[it.col()]) {
                    triplets.emplace_back(i, static_cast<int>(j), val * it.value());
                }
             }
         }

         return {rows, cols, triplets};
    }

    IntervalMatrix operator*(const Eigen::Matrix<zono_float, -1, -1>& A, const IntervalMatrix& B) {
        const Eigen::SparseMatrix<zono_float, Eigen::RowMajor> A_sp = A.sparseView();
        return A_sp * B;
    }

    IntervalMatrix IntervalMatrix::operator*(const IntervalMatrix& other) const
    {
        // input handling
        if (other.rows_ != this->cols_)
            throw std::invalid_argument("IntervalMatrix multiplication with IntervalMatrix: inconsistent dimensions");

        // dimensions
        size_t rows = this->rows_;
        size_t cols = other.cols_;

        // loop through to generate triplets
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (size_t i = 0; i < rows; ++i)
        {
            for (const auto& [k, val_a] : this->mat_[i])
            {
                for (const auto& [j, val_b] : other.mat_[k])
                {
                    triplets.emplace_back(static_cast<int>(i), static_cast<int>(j), val_a * val_b);
                }
            }
        }
        return {rows, cols, triplets};
    }

    void IntervalMatrix::operator*=(const IntervalMatrix& other)
    {
        *this = *this * other;
    }

    IntervalMatrix IntervalMatrix::operator+(const IntervalMatrix& other) const
    {
        // input handling
        if (other.rows_ != this->rows_ || other.cols_ != this->cols_)
            throw std::invalid_argument("IntervalMatrix addition with IntervalMatrix: inconsistent dimensions");

        // generate triplets and let constructor take care of adding them
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < static_cast<int>(this->rows_); ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, static_cast<int>(j), val);
            }
            for (const auto& [j, val] : other.mat_[i])
            {
                triplets.emplace_back(i, static_cast<int>(j), val);
            }
        }
        return {this->rows_, this->cols_, triplets};
    }

    void IntervalMatrix::operator+=(const IntervalMatrix& other)
    {
        *this = *this + other;
    }

    IntervalMatrix IntervalMatrix::operator-(const IntervalMatrix& other) const
    {
        // input handling
        if (other.rows_ != this->rows_ || other.cols_ != this->cols_)
            throw std::invalid_argument("IntervalMatrix subtraction with IntervalMatrix: inconsistent dimensions");

        // generate triplets and let constructor take care of adding them
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < static_cast<int>(this->rows_); ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, static_cast<int>(j), val);
            }
            for (const auto& [j, val] : other.mat_[i])
            {
                triplets.emplace_back(i, static_cast<int>(j), val * (-one));
            }
        }
        return {this->rows_, this->cols_, triplets};
    }

    IntervalMatrix IntervalMatrix::operator+(const Interval& interval) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val += interval;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator+=(const Interval& interval)
    {
        *this = *this + interval;
    }

    IntervalMatrix IntervalMatrix::operator+(zono_float alpha) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val += alpha;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator+=(zono_float alpha)
    {
        *this = *this + alpha;
    }

    IntervalMatrix operator+(zono_float alpha, const IntervalMatrix& A) {
         return A + alpha;
    }

    IntervalMatrix operator-(const Interval& interval, const IntervalMatrix& mat)
    {
        return -mat + interval;
    }

    void IntervalMatrix::operator-=(const IntervalMatrix& other)
    {
        *this = *this - other;
    }

    IntervalMatrix IntervalMatrix::operator-(const Interval& interval) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val -= interval;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator-=(const Interval& interval)
    {
        *this = *this - interval;
    }

    IntervalMatrix IntervalMatrix::operator-(zono_float alpha) const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val -= alpha;
            }
        }
        return mat;
    }

    void IntervalMatrix::operator-=(zono_float alpha)
    {
        *this = *this - alpha;
    }

    IntervalMatrix operator-(zono_float alpha, const IntervalMatrix& A)
    {
        return -A + alpha;
    }

    IntervalMatrix operator/(zono_float alpha, const IntervalMatrix& A)
    {
        IntervalMatrix mat(A);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val = alpha / val;
            }
        }
        return mat;
    }

    IntervalMatrix operator/(const Interval& interval, const IntervalMatrix& A)
    {
        IntervalMatrix mat(A);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val = interval / val;
            }
        }
        return mat;
    }

    IntervalMatrix IntervalMatrix::operator-() const
    {
        IntervalMatrix mat(*this);
        for (int i=0; i<static_cast<int>(mat.rows()); ++i)
        {
            for (auto& [k, val] : mat.mat_[i])
            {
                val = -val;
            }
        }
        return mat;
    }

    Eigen::SparseMatrix<zono_float> IntervalMatrix::center() const
    {
        std::vector<Eigen::Triplet<zono_float>> triplets;
        for (int k = 0; k < static_cast<int>(this->rows_); ++k)
        {
            for (const auto& [col, val] : this->mat_[k])
            {
                triplets.emplace_back(k, static_cast<int>(col), val.center());
            }
        }
        Eigen::SparseMatrix<zono_float> center_mat(static_cast<Eigen::Index>(this->rows_),
                                                   static_cast<Eigen::Index>(this->cols_));
        center_mat.setFromTriplets(triplets.begin(), triplets.end());
        return center_mat;
    }

    Eigen::SparseMatrix<zono_float> IntervalMatrix::diam() const
    {
        std::vector<Eigen::Triplet<zono_float>> triplets;
        for (int k = 0; k < static_cast<int>(this->rows_); ++k)
        {
            for (const auto& [col, val] : this->mat_[k])
            {
                triplets.emplace_back(k, static_cast<int>(col), val.width());
            }
        }
        Eigen::SparseMatrix<zono_float> diam_mat(static_cast<Eigen::Index>(this->rows_),
                                                 static_cast<Eigen::Index>(this->cols_));
        diam_mat.setFromTriplets(triplets.begin(), triplets.end());
        return diam_mat;
    }

    IntervalMatrix IntervalMatrix::radius() const
    {
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int i = 0; i < static_cast<int>(this->rows_); ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, static_cast<int>(j), val.radius());
            }
        }
        return {this->rows_, this->cols_, triplets};
    }

    zono_float IntervalMatrix::width() const
    {
        zono_float w = 0;
        for (size_t i = 0; i < this->rows_; i++)
        {
            for (const auto& [col, val] : this->mat_[i])
            {
                w = std::max(w, val.width());
            }
        }
        return w;
    }

    std::string IntervalMatrix::print() const
    {
        std::stringstream ss;
        ss << "IntervalMatrix: " << std::endl;
        ss << " rows = " << this->rows_ << std::endl;
        ss << " cols = " << this->cols_ << std::endl;
        ss << " elements: " << std::endl;
        for (size_t i = 0; i < this->rows_; i++)
        {
            for (const auto& [col, val] : this->mat_[i])
            {
                ss << "  (" << i << ", " << col << ", " << val << ")" << std::endl;
            }
        }
        return ss.str();
    }

    std::ostream& operator<<(std::ostream& os, const IntervalMatrix& interval_matrix)
    {
        os << interval_matrix.print();
        return os;
    }

}