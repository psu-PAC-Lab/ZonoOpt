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

        // in-place operators (type-preserving overrides)
        using HybZono::operator+=;  // restore vector overload hidden by declarations below
        using HybZono::operator-=;  // restore vector overload hidden by declarations below
        using HybZono::operator*=;  // restore scalar overload hidden by declarations below

        /**
         * @brief In-place Minkowski sum with another Point (translation).
         * Use operator+ for Zono/ConZono/HybZono/Box arguments, which return a wider type.
         */
        void operator+=(Point& other);

        /** @brief Deleted: use operator+ instead — result type would be Zono or wider. */
        void operator+=(HybZono& other) = delete;

        /** @brief Deleted: use operator+ instead — result type would be Zono. */
        void operator+=(const Box& box) = delete;

        /** @brief Deleted: use operator- instead — Point ⊖ Zono does not yield a Point. */
        void operator-=(Zono& other) = delete;

        /** @brief Deleted: use operator- instead — Point ⊖ Box does not yield a Point. */
        void operator-=(const Box& box) = delete;

        /**
         * @brief In-place Cartesian product with another Point.
         * Use operator* for Zono/ConZono/HybZono/Box arguments, which return a wider type.
         */
        void operator*=(Point& other);

        /** @brief Deleted: use operator* instead — result type would be Zono or wider. */
        void operator*=(HybZono& other) = delete;

        /** @brief Deleted: use operator* instead — result type would be Zono. */
        void operator*=(const Box& box) = delete;

        // do nothing methods
        std::unique_ptr<HybZono> remove_redundancy(int) const override { return std::unique_ptr<HybZono>(this->clone()); }

        void convert_form() override
        {
            /* do nothing */
        }

        std::unique_ptr<ConZono> constraint_reduction() const override
        {
            return std::make_unique<Point>(*this);
        }

    protected:
        Eigen::Vector<zono_float, -1> do_optimize_over(
            const Eigen::SparseMatrix<zono_float>&, const Eigen::Vector<zono_float, -1>&, zono_float,
            const SolverSettings&, std::shared_ptr<OptSolution>* sol,
            const WarmStartParams&) const override;

        Eigen::Vector<zono_float, -1> do_project_point(const Eigen::Vector<zono_float, -1>& x,
                                                       const SolverSettings&, std::shared_ptr<OptSolution>* sol,
                                                       const WarmStartParams&) const override;

        zono_float do_support(const Eigen::Vector<zono_float, -1>& d, const SolverSettings&,
                              std::shared_ptr<OptSolution>* sol, const WarmStartParams&) override;

        bool do_contains_point(const Eigen::Vector<zono_float, -1>& x, const SolverSettings&,
                               std::shared_ptr<OptSolution>* sol, const WarmStartParams&) const override;

        Box do_bounding_box(const SolverSettings&, std::shared_ptr<OptSolution>* sol, const WarmStartParams&) override;
    
    private:
        void make_default_solution(std::shared_ptr<OptSolution>* sol) const;
    
    };
} // namespace ZonoOpt

#endif
