#ifndef ZONOOPT_ZONO_HPP_
#define ZONOOPT_ZONO_HPP_

/**
 * @file Zono.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Zonotope class for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "ConZono.hpp"

namespace ZonoOpt
{
    using namespace detail;

    /**
     * @brief Zonotope class.
     *
     * A zonotope is defined as:
     * Z = {G * xi + c | xi in [-1, 1]^nG}.
     * Equivalently, the following shorthand can be used: Z = <G, c>.
     * Optionally, in 0-1 form, the factors are xi in [0, 1]^nG.
     * The set dimension is n, and the number of equality constraints is nC.
     *
     */
    class Zono : public ConZono
    {
    public:
        // constructors
        /**
         * @brief Default constructor for Zono class
         *
         */
        Zono() { sharp = true; }

        /**
         * @brief Zono constructor
         *
         * @param G generator matrix
         * @param c center
         * @param zero_one_form true if set is in 0-1 form
         */
        Zono(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
             const bool zero_one_form = false);

        // virtual destructor
        ~Zono() override = default;

        // set method
        /**
         * @brief Reset zonotope object with the given parameters.
         * 
         * @param G generator matrix
         * @param c center
         * @param zero_one_form true if set is in 0-1 form
         */
        void set(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                 bool zero_one_form = false);

        /**
         * @brief Clone method for polymorphic behavior.
         */
        HybZono* clone() const override;

        /**
         * @brief Perform zonotope order reduction
         *
         * @param n_o desired order, must be greater than or equal to the dimension of the set
         * @return zonotope with order n_o
         */
        std::unique_ptr<Zono> reduce_order(int n_o);

        /**
         * @brief Get volume of zonotope
         * @return volume
         *
         * Reference: Gover and Krikorian 2010, "Determinants and the volumes of parallelotopes and zonotopes"
         * Requires nG choose n determinant computations.
         */
        zono_float get_volume();

        /**
         * @brief get center of zonotope
         * @return center
         */
        Eigen::Vector<zono_float, -1> get_center();

        // generator conversion between [-1,1] and [0,1]
        void convert_form() override;

        // display methods
        std::string print() const override;

    protected:
        bool do_is_empty(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&) const override;

        Box do_bounding_box(const OptSettings&, std::shared_ptr<OptSolution>*, const WarmStartParams&) override;

        zono_float do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings&,
                              std::shared_ptr<OptSolution>*, const WarmStartParams&) override;
    };

    // forward declarations
    /**
     * @brief Builds a zonotope from a Box object.
     *
     * @param box Box object (vector of intervals)
     * @return zonotope
     * @ingroup ZonoOpt_SetupFunctions
     */
    std::unique_ptr<Zono> interval_2_zono(const Box& box);

    /**
     * @brief Builds a 2D regular zonotope with a given radius and number of sides.
     *
     * @param radius radius of the zonotope
     * @param n_sides number of sides (must be an even number >= 4)
     * @param outer_approx flag to do an outer approximation instead of an inner approximation
     * @param c center vector
     * @return zonotope
     * @ingroup ZonoOpt_SetupFunctions
     */
    std::unique_ptr<Zono> make_regular_zono_2D(zono_float radius, int n_sides, bool outer_approx = false,
                                               const Eigen::Vector<zono_float, 2>& c = Eigen::Vector<
                                                   zono_float, 2>::Zero());
} // namespace ZonoOpt

#endif
