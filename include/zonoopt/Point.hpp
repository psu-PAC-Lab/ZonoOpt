#ifndef ZONOOPT_POINT_HPP_
#define ZONOOPT_POINT_HPP_

/**
 * @file Point.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Point class for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "Zono.hpp"

namespace ZonoOpt
{
    using namespace detail;

    /**
     * @brief Point class.
     *
     * A point is defined entirely by the center vector c.
     */
    class Point final : public Zono
    {
    public:
        // constructor
        /**
         * @brief Default constructor for Point class
         *
         */
        Point() { sharp = true; }

        /**
         * @brief Point constructor
         *
         * @param c center
         */
        explicit Point(const Eigen::Vector<zono_float, -1>& c);

        // set method
        /**
         * @brief Reset point object with the given parameters.
         * 
         * @param c center
         */
        void set(const Eigen::Vector<zono_float, -1>& c);

        HybZono* clone() const override;

        // display methods
        std::string print() const override;

        // do nothing methods
        bool remove_redundancy(int) override { return false; }

        void convert_form() override
        {
            /* do nothing */
        }

    protected:
        Eigen::Vector<zono_float, -1> do_optimize_over(
            const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, zono_float,
            const OptSettings&, std::shared_ptr<OptSolution>*,
            const WarmStartParams&) const override;

        Eigen::Vector<zono_float, -1> do_project_point(const Eigen::Vector<zono_float, -1>& x,
                                                       const OptSettings&, std::shared_ptr<OptSolution>*,
                                                       const WarmStartParams&) const override;

        zono_float do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings&,
                              std::shared_ptr<OptSolution>*, const WarmStartParams&) override;

        bool do_contains_point(const Eigen::Vector<zono_float, -1>& x, const OptSettings&,
                               std::shared_ptr<OptSolution>*, const WarmStartParams&) const override;

        Box do_bounding_box(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&) override;
    };
} // namespace ZonoOpt

#endif
