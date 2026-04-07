#include "ZonoOpt.hpp"

namespace ZonoOpt {
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
            this->x_lb(i) = vals[i].lower();
            this->x_ub(i) = vals[i].upper();
        }
    }

    Box::Box(const Eigen::Vector<Interval, -1>& vals)
    {
        x_lb.resize(vals.size());
        x_ub.resize(vals.size());

        for (Eigen::Index i = 0; i < vals.size(); i++)
        {
            this->x_lb(i) = vals(i).lower();
            this->x_ub(i) = vals(i).upper();
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

    std::vector<Interval> Box::to_array() const
    {
        std::vector<Interval> vals;
        for (int i = 0; i < x_lb.size(); ++i)
        {
            vals.emplace_back(x_lb(i), x_ub(i));
        }
        return vals;
    }

    Interval Box::get_element(const int i) const
    {
        if (i >= x_lb.size())
            throw std::out_of_range("Index out of range");
        return {x_lb(i), x_ub(i)};
    }

    void Box::set_element(const int i, const Interval& val)
    {
        this->x_lb(i) = val.lower();
        this->x_ub(i) = val.upper();
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
        for (int i = 0; i < x_lb.size(); i++)
        {
            w = std::max(w, this->get_element(i).width());
        }
        return w;
    }

    Eigen::Vector<zono_float, -1> Box::center() const
    {
        Eigen::Vector<zono_float, -1> c(this->size());
        for (int i = 0; i < static_cast<int>(this->size()); i++)
        {
            c(i) = get_element(i).center();
        }
        return c;
    }

    Box Box::radius() const
    {
        Box out(this->size());
        for (int i = 0; i < static_cast<int>(this->size()); i++)
        {
            out.set_element(i, get_element(i).radius());
        }
        return out;
    }

    Box Box::intersect(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box intersection: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i).intersect(other.get_element(i)));
        }
        return out;
    }

    Box Box::interval_hull(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box interval hull: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i).interval_hull(other.get_element(i)));
        }
        return out;
    }

    bool Box::contains(const Eigen::Vector<zono_float, -1>& v)
    {
        if (v.size() != static_cast<Eigen::Index>(this->size()))
            throw std::invalid_argument("Box contains: inconsistent dimensions");
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            if (!get_element(i).contains(v(i)))
                return false;
        }
        return true;
    }

    bool Box::contains_set(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box contains_set: inconsistent dimensions");
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            if (!get_element(i).contains_set(other.get_element(i)))
                return false;
        }
        return true;
    }

    Box Box::operator+(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box addition: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i) + other.get_element(i));
        }
        return out;
    }

    void Box::operator+=(const Box& other)
    {
        *this = *this + other;
    }

    Box Box::operator+(const Eigen::Vector<zono_float, -1>& v) const
    {
        if (v.size() != static_cast<Eigen::Index>(this->size()))
            throw std::invalid_argument("Box addition with vector: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i) + v(i));
        }
        return out;
    }

    void Box::operator+=(const Eigen::Vector<zono_float, -1>& v)
    {
        *this = *this + v;
    }

    Box operator+(const Eigen::Vector<zono_float, -1>& v, const Box& box)
    {
        return box + v;
    }

    Box Box::operator-(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box subtraction: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) - other.get_element(i));
        }
        return out;
    }

    void Box::operator-=(const Box& other)
    {
        *this = *this - other;
    }

    Box Box::operator-(const Eigen::Vector<zono_float, -1>& v) const
    {
        return *this + (-v);
    }

    void Box::operator-=(const Eigen::Vector<zono_float, -1>& v)
    {
        *this = *this - v;
    }

    Box operator-(const Eigen::Vector<zono_float, -1>& v, const Box& box)
    {
        return v + (-box);
    }

    void Box::operator*=(const Box& other)
    {
        *this = *this * other;
    }

    void Box::operator*=(zono_float alpha)
    {
        *this = *this * alpha;
    }

    Box operator*(zono_float alpha, const Box& box)
    {
        return box*alpha;
    }

    Box Box::operator*(const Eigen::Vector<zono_float, -1>& v) const
    {
        if (v.size() != static_cast<Eigen::Index>(this->size()))
            throw std::invalid_argument("Box multiplication with vector: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i) * v(i));
        }
        return out;
    }

    void Box::operator*=(const Eigen::Vector<zono_float, -1>& v)
    {
        *this = *this * v;
    }

    Box operator*(const Eigen::Vector<zono_float, -1>& v, const Box& box)
    {
        return box*v;
    }

    Box Box::operator*(const Interval& interval) const
    {
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, get_element(i) * interval);
        }
        return out;
    }

    Box operator*(const Interval& interval, const Box& box)
    {
        return box*interval;
    }

    void Box::operator*=(const Interval& interval)
    {
        *this = *this * interval;
    }

    Box operator*(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A, const Box& box)
    {
        return box.linear_map(A);
    }

    Box operator*(const Eigen::Matrix<zono_float, -1, -1>& A, const Box& box)
    {
        return box.linear_map(A);
    }

    Box Box::operator*(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box multiplication: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) * other.get_element(i));
        }
        return out;
    }

    Box Box::operator*(zono_float alpha) const
    {
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) * alpha);
        }
        return out;
    }

    Box Box::operator/(const Box& other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Box division: inconsistent dimensions");
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) / other.get_element(i));
        }
        return out;
    }

    void Box::operator/=(const Box& other)
    {
        *this = *this / other;
    }

    Box Box::operator/(zono_float alpha) const
    {
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) / alpha);
        }
        return out;
    }

    void Box::operator/=(zono_float alpha)
    {
        *this = *this / alpha;
    }

    Box operator/(zono_float alpha, const Box& box)
    {
        Box out = box;
        for (int i = 0; i < static_cast<int>(box.size()); ++i)
        {
            out.set_element(i, alpha / box.get_element(i));
        }
        return out;
    }

    Box Box::operator/(const Interval& interval) const
    {
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, this->get_element(i) / interval);
        }
        return out;
    }

    void Box::operator/=(const Interval& interval)
    {
        *this = *this / interval;
    }

    Box operator/(const Interval& interval, const Box& box)
    {
        Box out = box;
        for (int i = 0; i < static_cast<int>(box.size()); ++i)
        {
            out.set_element(i, interval / box.get_element(i));
        }
        return out;
    }

    Box Box::operator-() const
    {
        Box out = *this;
        for (int i = 0; i < static_cast<int>(this->size()); ++i)
        {
            out.set_element(i, -this->get_element(i));
        }
        return out;
    }

    Box Box::operator&(const Box& other) const
    {
        return this->intersect(other);
    }

    Box Box::operator|(const Box& other) const
    {
        return this->interval_hull(other);
    }

    bool Box::operator<=(const Box& other) const
    {
        return other.contains_set(*this);
    }

    bool Box::operator>=(const Box& other) const
    {
        return this->contains_set(other);
    }

    bool Box::operator==(const Box& other) const
    {
        return this->contains_set(other) && other.contains_set(*this);
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
            y.set_element(i, Interval(0, 0));
            for (int j = 0; j < A.cols(); j++)
            {
                y.set_element(i, y.get_element(i) + this->get_element(j)*A(i, j));
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
            y.set_element(i, Interval(0, 0));
            for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it(A, i); it; ++it)
            {
                y.set_element(i, y.get_element(i) + this->get_element(static_cast<int>(it.col()))*it.value());
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
            y = y + get_element(i)*x(i);
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
        for (int i = 0; i < static_cast<int>(x_lb.size()); i++)
        {
            ss << "  " << this->get_element(i) << std::endl;
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
                y = y + this->get_element(static_cast<int>(it.col()))*it.value();
            }

            // check validity
            if (!y.contains(b(k)))
                return false;
        }

        return true; // constraints valid
    }

    void Box::contract_backward(const Eigen::SparseMatrix<zono_float, Eigen::RowMajor>& A,
                                const Eigen::Vector<zono_float, -1>& b, const std::set<int>& constraints)
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

                for (Eigen::SparseMatrix<zono_float, Eigen::RowMajor>::InnerIterator it_inner(A, k); it_inner; ++it_inner)
                {
                    if (it_inner.col() != it_outer.col())
                    {
                        x = x + this->get_element(static_cast<int>(it_inner.col()))*(-it_inner.value());
                    }
                }

                // update interval
                this->set_element(static_cast<int>(it_outer.col()), this->get_element(static_cast<int>(it_outer.col())).intersect(x * (one/a_col)));
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

    MI_Box::MI_Box(const std::vector<Interval>& intervals, const std::pair<int, int>& idx_b, bool zero_one_form) : Box(intervals), idx_b(idx_b)
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
                if (this->get_element(ib).is_single_valued())
                    continue; // already fixed

                // test if low value is valid
                Box low_box(this->x_lb, this->x_ub);
                low_box.set_element(ib, Interval(this->bin_low, this->bin_low));
                low_box.contract_backward(A, b, constraints);
                const bool low_valid = low_box.contract_forward(A, b, constraints);

                // test if high value is valid
                Box high_box(this->x_lb, this->x_ub);
                high_box.set_element(ib, Interval(this->bin_high, this->bin_high));
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
}