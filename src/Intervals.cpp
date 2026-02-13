#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    Box::Box(const size_t size)
    {
        x_lb.resize(static_cast<Eigen::Index>(size));
        x_ub.resize(static_cast<Eigen::Index>(size));
    }

    Box::Box(const std::vector<Interval>& vals)
    {
        x_lb.resize(static_cast<Eigen::Index>(vals.size()));
        x_ub.resize(static_cast<Eigen::Index>(vals.size()));

        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(vals.size()); i++)
        {
            this->x_lb(i) = vals[i].y_min();
            this->x_ub(i) = vals[i].y_max();
        }
    }

    Box::Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub)
    {
        if (x_lb.size() != x_ub.size())
            throw std::invalid_argument("x_l and x_u must have the same size");
        this->x_lb = x_lb;
        this->x_ub = x_ub;
    }

    Box& Box::operator=(const Box& other)
    {
        if (this != &other)
        {
            this->x_lb = other.x_lb;
            this->x_ub = other.x_ub;
        }
        return *this;
    }

    Box::Box(const Box& other)
    {
        this->x_lb = other.x_lb;
        this->x_ub = other.x_ub;
    }

    IntervalView Box::operator[](const size_t i)
    {
        if (i >= static_cast<size_t>(x_lb.size()))
            throw std::out_of_range("Index out of range");
        return IntervalView(&x_lb(static_cast<Eigen::Index>(i)), &x_ub(static_cast<Eigen::Index>(i)));
    }

    Interval Box::operator[](const size_t i) const
    {
        if (i >= static_cast<size_t>(x_lb.size()))
            throw std::out_of_range("Index out of range");
        return Interval(x_lb(static_cast<Eigen::Index>(i)), x_ub(static_cast<Eigen::Index>(i)));
    }

    size_t Box::size() const
    {
        return x_lb.size();
    }

    void Box::project(Eigen::Ref<Eigen::Vector<zono_float, -1>> x) const
    {
        if (x.size() != x_lb.size())
            throw std::invalid_argument("x must have the same size as the Box");
        x = x.cwiseMax(x_lb).cwiseMin(x_ub);
    }

    Box* Box::clone() const
    {
        return new Box(*this);
    }

    zono_float Box::width() const
    {
        zono_float w = 0;
        for (Eigen::Index i = 0; i < x_lb.size(); i++)
        {
            w = std::max(w, (*this)[i].width());
        }
        return w;
    }

    Eigen::Vector<zono_float, -1> Box::center() const
    {
        Eigen::Vector<zono_float, -1> c(this->size());
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(this->size()); i++)
        {
            c(i) = (*this)[i].center();
        }
        return c;
    }

    Box Box::radius() const
    {
        Box out(this->size());
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(this->size()); i++)
        {
            out[i] = (*this)[i].radius();
        }
        return out;
    }

    Box Box::operator+(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box addition: inconsistent dimensions");
        Box out = *this;
        for (size_t i = 0; i < this->size(); ++i)
        {
            out[i] = (*this)[i] + other[i];
        }
        return out;
    }

    Box Box::operator-(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box subtraction: inconsistent dimensions");
        Box out = *this;
        for (size_t i = 0; i < this->size(); ++i)
        {
            out[i] = (*this)[i] - other[i];
        }
        return out;
    }

    Box Box::operator*(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box multiplication: inconsistent dimensions");
        Box out = *this;
        for (size_t i = 0; i < this->size(); ++i)
        {
            out[i] = (*this)[i] * other[i];
        }
        return out;
    }

    Box Box::operator*(zono_float alpha) const
    {
        Box out = *this;
        for (size_t i = 0; i < this->size(); ++i)
        {
            out[i] = (*this)[i] * alpha;
        }
        return out;
    }

    Box Box::operator/(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box division: inconsistent dimensions");
        Box out = *this;
        for (size_t i = 0; i < this->size(); ++i)
        {
            out[i] = (*this)[i] / other[i];
        }
        return out;
    }

    bool Box::contract(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                       const Eigen::Vector<zono_float, -1>& b, const int iter)
    {
        if (iter <= 0)
            throw std::invalid_argument("iter must be positive");

        // contract over all constraints
        std::set<int> constraints;
        for (int i = 0; i < A.rows(); i++)
        {
            constraints.insert(i);
        }

        // run contractor
        return contract_helper(A, b, iter, constraints);
    }

    bool Box::contract_subset(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                              const Eigen::Vector<zono_float, -1>& b, int iter,
                              const Eigen::SparseMatrix<zono_float>& A, const std::set<int>& inds,
                              int tree_search_depth)
    {
        if (iter <= 0)
            throw std::invalid_argument("iter must be positive");
        if (A_rm.rows() != A.rows() || A_rm.cols() != A.cols())
            throw std::invalid_argument("A_rm must equal A");

        // get affected constraints
        std::set<int> all_constraints, all_vars;

        get_vars_cons(A, A_rm, all_constraints, all_vars, inds, 0, tree_search_depth);

        // run contractor
        return contract_helper(A_rm, b, iter, all_constraints);
    }

    Box Box::linear_map(const Eigen::Matrix<zono_float, -1, -1>& A) const
    {
        // input handling
        if (A.cols() != x_lb.size())
            throw std::invalid_argument("Matrix A must have the same number of columns as the size of the Box");

        // declare
        Box y(A.rows());

        // linear map
        for (int i = 0; i < A.rows(); i++)
        {
            y[i] = Interval(0, 0);
            for (int j = 0; j < A.cols(); j++)
            {
                y[i].add_assign(y[i], ((*this)[j] * A(i, j)).as_view());
            }
        }
        return y;
    }

    Box Box::linear_map(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A) const
    {
        // input handling
        if (A.cols() != x_lb.size())
            throw std::invalid_argument("Matrix A must have the same number of columns as the size of the Box");

        // declare
        Box y(A.rows());

        // linear map
        for (int i = 0; i < A.rows(); i++)
        {
            y[i] = Interval(0, 0);
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it)
            {
                y[i].add_assign(y[i], ((*this)[it.col()] * it.value()).as_view());
            }
        }
        return y;
    }

    Interval Box::dot(const Eigen::Vector<zono_float, -1>& x) const
    {
        // input handling
        if (x.size() != x_lb.size())
            throw std::invalid_argument("Vector x must have the same size as the Box");

        // declare
        Interval y(0, 0);

        // linear map
        for (int i = 0; i < this->x_lb.size(); i++)
            y.add_assign(y, ((*this)[i] * x(i)));
        return y;
    }

    void Box::permute(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P)
    {
        this->x_lb = P * this->x_lb;
        this->x_ub = P * this->x_ub;
    }

    std::string Box::print() const
    {
        std::stringstream ss;
        ss << "Box: " << std::endl;
        for (Eigen::Index i = 0; i < x_lb.size(); i++)
        {
            ss << "  " << (*this)[i] << std::endl;
        }
        return ss.str();
    }

    std::ostream& operator<<(std::ostream& os, const Box& box)
    {
        os << box.print();
        return os;
    }

    bool Box::contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                              const Eigen::Vector<zono_float, -1>& b, const int iter,
                              const std::set<int>& constraints)
    {
        for (int i = 0; i < iter; i++)
        {
            // forward
            if (!contract_forward(A, b, constraints))
                return false;

            // backward
            contract_backward(A, b, constraints);
        }

        return true; // constraints valid
    }

    bool Box::contract_forward(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                               const Eigen::Vector<zono_float, -1>& b, const std::set<int>& constraints)
    {
        // loop through constraints
        for (const int k : constraints)
        {
            // forward propagate
            Interval y(0, 0);
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
            {
                y.add_assign(y, (*this)[it.col()].to_interval() * it.value());
            }

            // check validity
            if (!y.contains(b(k)))
                return false;
        }

        return true; // constraints valid
    }

    void Box::contract_backward(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                const Eigen::Vector<double, -1>& b, const std::set<int>& constraints)
    {
        // loop through constraints
        for (const int k : constraints)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_outer(A, k); it_outer; ++it_outer)
            {
                auto x = Interval(b(k), b(k)); // init
                const zono_float a_col = it_outer.value();
                if (std::abs(a_col) < zono_eps)
                    continue; // skip zero coefficient

                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_inner(A, k); it_inner; ++
                     it_inner)
                {
                    if (it_inner.col() != it_outer.col())
                    {
                        x.add_assign(x, (*this)[it_inner.col()].to_interval() * (-it_inner.value()));
                    }
                }

                // update interval
                (*this)[it_outer.col()].intersect_assign((*this)[it_outer.col()], (x * (one / a_col)).as_view());
            }
        }
    }

    void Box::get_vars_cons(const Eigen::SparseMatrix<zono_float>& A,
                            const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A_rm,
                            std::set<int>& constraints, std::set<int>& vars, const std::set<int>& new_vars, int depth,
                            int max_depth)
    {
        // immediately copy over constraints and vars
        vars.insert(new_vars.begin(), new_vars.end());

        // find new constraints
        std::set<int> new_constraints;
        for (const int i : new_vars)
        {
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it(A, i); it; ++it)
            {
                if (!constraints.count(static_cast<int>(it.row())))
                {
                    new_constraints.insert(static_cast<int>(it.row()));
                }
            }
        }

        // find new vars
        std::set<int> new_new_vars;
        for (const int i : new_constraints)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A_rm, i); it; ++it)
            {
                if (!vars.count(static_cast<int>(it.col())))
                {
                    new_new_vars.insert(static_cast<int>(it.col()));
                }
            }
        }

        // add new constraints to set
        constraints.insert(new_constraints.begin(), new_constraints.end());

        // recurse if able
        depth++;
        if (depth < max_depth && new_new_vars.empty())
            get_vars_cons(A, A_rm, constraints, vars, new_new_vars, depth, max_depth);
    }

    MI_Box::MI_Box(const Eigen::Vector<zono_float, -1>& x_lb, const Eigen::Vector<zono_float, -1>& x_ub,
                   const std::pair<int, int>& idx_b, const bool zero_one_form) : Box(x_lb, x_ub), idx_b(idx_b)
    {
        // check that binary variables are in range
        if (idx_b.first < 0 || idx_b.first + idx_b.second > static_cast<int>(this->size()))
            throw std::out_of_range("Binary variable index out of range");

        // zero-one form handling
        this->bin_low = (zero_one_form) ? zero : -one;
        this->bin_high = one;
    }

    Box* MI_Box::clone() const
    {
        return new MI_Box(*this);
    }

    void MI_Box::project(Eigen::Ref<Eigen::Vector<zono_float, -1>> x) const
    {
        if (x.size() != x_lb.size())
            throw std::invalid_argument("x must have the same size as the Box");

        // first project onto box
        Box::project(x);

        // binary variable projection
        if (this->idx_b.second <= 0)
            return;

        const auto x_b = x.segment(idx_b.first, idx_b.second).array();
        const auto x_l_b = this->x_lb.segment(idx_b.first, idx_b.second).array();
        const auto x_u_b = this->x_ub.segment(idx_b.first, idx_b.second).array();

        const auto condition = x_b < (x_l_b + x_u_b) * 0.5;

        x.segment(this->idx_b.first, this->idx_b.second) = condition.select(x_l_b, x_u_b);
    }

    bool MI_Box::contract_helper(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                                 const Eigen::Vector<zono_float, -1>& b, const int iter,
                                 const std::set<int>& constraints)
    {
        // find participating binary variables
        std::set<int> bin_vars;
        for (const int k : constraints)
        {
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it)
            {
                if (it.col() >= this->idx_b.first && it.col() < this->idx_b.first + this->idx_b.second)
                {
                    bin_vars.insert(static_cast<int>(it.col()));
                }
            }
        }

        // contract with binary variable probing
        for (int i = 0; i < iter; ++i)
        {
            // binary variable probing
            for (const int ib : bin_vars)
            {
                if ((*this)[ib].is_single_valued())
                    continue; // already fixed

                // test if low value is valid
                Box low_box(this->x_lb, this->x_ub);
                low_box[ib] = Interval(this->bin_low, this->bin_low);
                low_box.contract_backward(A, b, constraints);
                const bool low_valid = low_box.contract_forward(A, b, constraints);

                // test if high value is valid
                Box high_box(this->x_lb, this->x_ub);
                high_box[ib] = Interval(this->bin_high, this->bin_high);
                high_box.contract_backward(A, b, constraints);
                const bool high_valid = high_box.contract_forward(A, b, constraints);

                // update interval based on validity
                if (low_valid && !high_valid)
                {
                    this->x_lb = low_box.x_lb;
                    this->x_ub = low_box.x_ub;
                }
                else if (high_valid && !low_valid)
                {
                    this->x_lb = high_box.x_lb;
                    this->x_ub = high_box.x_ub;
                }
                else if (!low_valid && !high_valid)
                {
                    return false; // both invalid
                }
            }
        }

        return true; // constraints valid
    }

    // interval matrices

    IntervalMatrix::IntervalMatrix(size_t rows, size_t cols, const std::vector<Eigen::Triplet<Interval>>& triplets)
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

    IntervalMatrix::IntervalMatrix(const Eigen::SparseMatrix<zono_float>& mat_lb,
                                   const Eigen::SparseMatrix<zono_float>& mat_ub)
    {
        if (mat_lb.rows() != mat_ub.rows() || mat_lb.cols() != mat_ub.cols())
            throw std::invalid_argument("IntervalMatrix: mat_lb and mat_ub must have the same dimensions");

        // get triplets
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (int k = 0; k < mat_lb.outerSize(); k++)
        {
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it_lb(mat_lb, k); it_lb; ++it_lb)
            {
                triplets.emplace_back(it_lb.row(), it_lb.col(), Interval(it_lb.value(), zero));
            }
        }
        for (int k = 0; k < mat_ub.outerSize(); k++)
        {
            for (Eigen::SparseMatrix<zono_float>::InnerIterator it_ub(mat_ub, k); it_ub; ++it_ub)
            {
                triplets.emplace_back(it_ub.row(), it_ub.col(), Interval(zero, it_ub.value()));
            }
        }

        // let logic in triplet constructor take care of combining
        *this = IntervalMatrix(static_cast<size_t>(mat_lb.rows()), static_cast<size_t>(mat_lb.cols()), triplets);
    }

    Box IntervalMatrix::operator*(const Eigen::Vector<zono_float, -1>& v) const
    {
        // input handling
        if (v.size() != static_cast<Eigen::Index>(this->cols_))
            throw std::invalid_argument("IntervalMatrix multiplication with box: inconsistent dimensions");

        // declare
        Box y(this->rows_);

        // linear map
        for (size_t i = 0; i < this->rows_; i++)
        {
            y[i] = Interval(0, 0);
            for (const auto& [j, val] : this->mat_[i])
            {
                y[i].add_assign(y[i], (val * v[static_cast<Eigen::Index>(j)]).as_view());
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
        for (size_t i = 0; i < this->rows_; i++)
        {
            y[i] = Interval(0, 0);
            for (const auto& [j, val] : this->mat_[i])
            {
                y[i].add_assign(y[i], (box[j] * val).as_view());
            }
        }
        return y;
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
        for (size_t i = 0; i < rows; ++i)
        {
            for (const auto& [k, val] : this->mat_[i])
            {
                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_A(
                         A, static_cast<Eigen::Index>(k)); it_A; ++it_A)
                {
                    triplets.emplace_back(i, it_A.col(), val * it_A.value());
                }
            }
        }
        return {rows, cols, triplets};
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
                    triplets.emplace_back(i, j, val_a * val_b);
                }
            }
        }
        return {rows, cols, triplets};
    }

    IntervalMatrix IntervalMatrix::operator+(const IntervalMatrix& other) const
    {
        // input handling
        if (other.rows_ != this->rows_ || other.cols_ != this->cols_)
            throw std::invalid_argument("IntervalMatrix addition with IntervalMatrix: inconsistent dimensions");

        // generate triplets and let constructor take care of adding them
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, j, val);
            }
            for (const auto& [j, val] : other.mat_[i])
            {
                triplets.emplace_back(i, j, val);
            }
        }
        return {this->rows_, this->cols_, triplets};
    }

    IntervalMatrix IntervalMatrix::operator-(const IntervalMatrix& other) const
    {
        // input handling
        if (other.rows_ != this->rows_ || other.cols_ != this->cols_)
            throw std::invalid_argument("IntervalMatrix subtraction with IntervalMatrix: inconsistent dimensions");

        // generate triplets and let constructor take care of adding them
        std::vector<Eigen::Triplet<Interval>> triplets;
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, j, val);
            }
            for (const auto& [j, val] : other.mat_[i])
            {
                triplets.emplace_back(i, j, val * (-one));
            }
        }
        return {this->rows_, this->cols_, triplets};
    }

    Eigen::SparseMatrix<zono_float> IntervalMatrix::center() const
    {
        std::vector<Eigen::Triplet<zono_float>> triplets;
        for (size_t k = 0; k < this->rows_; ++k)
        {
            for (const auto& [col, val] : this->mat_[k])
            {
                triplets.emplace_back(k, col, val.center());
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
        for (size_t k = 0; k < this->rows_; ++k)
        {
            for (const auto& [col, val] : this->mat_[k])
            {
                triplets.emplace_back(k, col, val.width());
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
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (const auto& [j, val] : this->mat_[i])
            {
                triplets.emplace_back(i, j, val.radius());
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
