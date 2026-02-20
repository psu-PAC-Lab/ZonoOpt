#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    Point::Point(const Eigen::Vector<zono_float, -1>& c)
    {
        set(c);
        sharp = true;
    }

    HybZono* Point::clone() const
    {
        return new Point(*this);
    }

    void Point::set(const Eigen::Vector<zono_float, -1>& c)
    {
        // point parameters
        this->c = c;
        this->n = static_cast<int>(this->c.size());

        // hybzono parameters
        this->G.resize(this->n, 0);
        this->nG = 0;
        this->nGc = this->nG;
        this->nGb = 0;
        this->nC = 0;
        this->Gc = this->G;
        this->Gb.resize(this->n, 0);
        this->A.resize(0, this->nG);
        this->Ac = this->A;
        this->Ab.resize(0, 0);
        this->b.resize(0);
        this->zero_one_form = false;
    }

    std::string Point::print() const
    {
        std::stringstream ss;
        ss << "Point: " << std::endl;
        ss << "n: " << this->n << std::endl;
        ss << "c: " << this->c;
        return ss.str();
    }

    Eigen::Vector<zono_float, -1> Point::do_optimize_over(
        const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, zono_float,
        const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&) const
    {
        return this->c;
    }

    Eigen::Vector<zono_float, -1> Point::do_project_point(const Eigen::Vector<zono_float, -1>& x,
                                                          const OptSettings&, std::shared_ptr<OptSolution>*,
                                                          const WarmStartParams&) const
    {
        // check dimensions
        if (this->n != x.size())
        {
            throw std::invalid_argument("Point projection: inconsistent dimensions.");
        }

        return this->c;
    }

    zono_float Point::do_support(const Eigen::Vector<zono_float, -1>& d,
                                 const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&)
    {
        // check dimensions
        if (this->n != d.size())
        {
            throw std::invalid_argument("Support: inconsistent dimensions.");
        }

        return this->c.dot(d);
    }

    bool Point::do_contains_point(const Eigen::Vector<zono_float, -1>& x,
                                  const OptSettings&, std::shared_ptr<OptSolution>*,
                                  const WarmStartParams&) const
    {
        if (this->n != x.size())
            throw std::invalid_argument("Contains point: inconsistent dimensions");

        const zono_float dist = (x - this->c).norm();
        return dist < zono_eps;
    }

    Box Point::do_bounding_box(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&)
    {
        return {this->c, this->c};
    }
}
