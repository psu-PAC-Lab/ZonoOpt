#include "ZonoOpt.hpp"

namespace ZonoOpt
{
    using namespace detail;

    EmptySet::EmptySet(const int n)
    {
        this->n = n;

        // hybzono parameters
        this->c.resize(this->n);
        this->c.setConstant(std::numeric_limits<zono_float>::quiet_NaN());
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

    HybZono* EmptySet::clone() const
    {
        return new EmptySet(*this);
    }

    std::string EmptySet::print() const
    {
        std::stringstream ss;
        ss << "EmptySet: " << std::endl;
        ss << "  n: " << this->n;
        return ss.str();
    }

    Eigen::Vector<zono_float, -1> EmptySet::do_optimize_over(
        const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, zono_float,
        const OptSettings&, std::shared_ptr<OptSolution>* solution,
        const WarmStartParams&) const
    {
        if (solution)
        {
            (*solution)->infeasible = true;
        }
        return Eigen::Vector<zono_float, -1>::Constant(this->n, std::numeric_limits<zono_float>::quiet_NaN());
    }

    Eigen::Vector<zono_float, -1> EmptySet::do_project_point(const Eigen::Vector<zono_float, -1>&, const OptSettings&,
                                                             std::shared_ptr<OptSolution>* solution,
                                                             const WarmStartParams&) const
    {
        if (solution)
        {
            (*solution)->infeasible = true;
        }
        return Eigen::Vector<zono_float, -1>::Constant(this->n, std::numeric_limits<zono_float>::quiet_NaN());
    }

    zono_float EmptySet::do_support(const Eigen::Vector<zono_float, -1>&, const OptSettings&,
                                    std::shared_ptr<OptSolution>* solution,
                                    const WarmStartParams&)
    {
        if (solution)
        {
            (*solution)->infeasible = true;
        }
        return std::numeric_limits<zono_float>::quiet_NaN();
    }

    bool EmptySet::do_contains_point(const Eigen::Vector<zono_float, -1>&, const OptSettings&,
                                     std::shared_ptr<OptSolution>*,
                                     const WarmStartParams&) const
    {
        return false;
    }

    Box EmptySet::do_bounding_box(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&)
    {
        const Eigen::Vector<zono_float, -1> x_l = Eigen::Vector<zono_float, -1>::Constant(
            this->n, std::numeric_limits<zono_float>::infinity());
        const Eigen::Vector<zono_float, -1> x_u = -Eigen::Vector<zono_float, -1>::Constant(
            this->n, std::numeric_limits<zono_float>::infinity());
        return {x_l, x_u};
    }

    bool EmptySet::do_is_empty(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&) const
    {
        return true;
    }

    std::unique_ptr<HybZono> EmptySet::do_complement(zono_float delta_m, bool, const OptSettings&,
                                                     std::shared_ptr<OptSolution>*,
                                                     int, int)
    {
        const zono_float m = delta_m + 1; // box width
        const Eigen::Vector<zono_float, -1> x_l = -Eigen::Vector<zono_float, -1>::Constant(this->n, m);
        const Eigen::Vector<zono_float, -1> x_u = Eigen::Vector<zono_float, -1>::Constant(this->n, m);
        const Box box(x_l, x_u);
        return interval_2_zono(box);
    }
}
