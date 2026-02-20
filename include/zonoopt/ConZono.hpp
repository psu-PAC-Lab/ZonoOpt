#ifndef ZONOOPT_CONZONO_HPP_
#define ZONOOPT_CONZONO_HPP_

/**
 * @file ConZono.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Constrained zonotope class for ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include "HybZono.hpp"

namespace ZonoOpt
{
    using namespace detail;

    /**
     * @brief Constrained zonotope class.
     *
     * A constrained zonotope is defined as:
     * Z = {G * xi + c | A * xi = b, xi in [-1, 1]^nG}.
     * Equivalently, the following shorthand can be used: Z = <G, c, A, b>.
     * Optionally, in 0-1 form, the factors are xi in [0, 1]^nG.
     * The set dimension is n, and the number of equality constraints is nC.
     *
     */
    class ConZono : public HybZono
    {
    public:
        // constructors

        /**
         * @brief Default constructor for ConZono class
         *
         */
        ConZono() { sharp = true; }

        /**
         * @brief ConZono constructor
         *
         * @param G generator matrix
         * @param c center
         * @param A constraint matrix
         * @param b constraint vector
         * @param zero_one_form true if set is in 0-1 form
         */
        ConZono(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                const bool zero_one_form = false);

        // virtual destructor
        ~ConZono() override = default;

        /**
         * @brief Clone method for polymorphic behavior.
         */
        HybZono* clone() const override;

        // set method
        /**
         * @brief Reset constrained zonotope object with the given parameters.
         *
         * @param G generator matrix
         * @param c center
         * @param A constraint matrix
         * @param b constraint vector
         * @param zero_one_form true if set is in 0-1 form
         */
        void set(const Eigen::SparseMatrix<zono_float>& G, const Eigen::Vector<zono_float, -1>& c,
                 const Eigen::SparseMatrix<zono_float>& A, const Eigen::Vector<zono_float, -1>& b,
                 bool zero_one_form = false);

        /**
         * @brief Execute constraint reduction algorithm from Scott et. al. 2016
         *
         * Removes one constraint and one generator from the constrained zonotope.
         * The resulting set is an over-approximation of the original set.
         */
        virtual void constraint_reduction();

        // generator conversion between [-1,1] and [0,1]
        void convert_form() override;

        // over-approximate as zonotope

        /**
         * @brief Compute outer approximation of constrained zonotope as zonotope using SVD
         * @return Zonotope over-approximation
         */
        virtual std::unique_ptr<Zono> to_zono_approx() const;

        // display methods
        std::string print() const override;

    protected:
        OptSolution qp_opt(const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q,
                           zono_float c, const Eigen::SparseMatrix<zono_float>& A,
                           const Eigen::Vector<zono_float, -1>& b,
                           const OptSettings& settings = OptSettings(),
                           std::shared_ptr<OptSolution>* solution = nullptr,
                           const WarmStartParams& warm_start_params = WarmStartParams()) const;

        Eigen::Vector<zono_float, -1> do_optimize_over(
            const Eigen::SparseMatrix<zono_float>& P, const Eigen::Vector<zono_float, -1>& q, const zono_float c,
            const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
            const WarmStartParams&) const override;

        Eigen::Vector<zono_float, -1> do_project_point(const Eigen::Vector<zono_float, -1>& x,
                                                       const OptSettings& settings,
                                                       std::shared_ptr<OptSolution>* solution,
                                                       const WarmStartParams& warm_start_params) const override;

        bool do_is_empty(const OptSettings& settings, std::shared_ptr<OptSolution>* solution,
                         const WarmStartParams& warm_start_params) const override;

        zono_float do_support(const Eigen::Vector<zono_float, -1>& d, const OptSettings& settings,
                              std::shared_ptr<OptSolution>* solution,
                              const WarmStartParams& warm_start_params) override;

        bool do_contains_point(const Eigen::Vector<zono_float, -1>& x, const OptSettings& settings,
                               std::shared_ptr<OptSolution>* solution,
                               const WarmStartParams& warm_start_params) const override;

        Box do_bounding_box(const OptSettings& settings, std::shared_ptr<OptSolution>*,
                            const WarmStartParams& warm_start_params) override;

        std::unique_ptr<HybZono> do_complement(zono_float delta_m, bool, const OptSettings&,
                                               std::shared_ptr<OptSolution>*, int, int) override;
    };

    // forward declarations

    /**
    * @brief Builds a constrained zonotope from a vertex representation polytope.
    *
    * @param Vpoly vertices of V-rep polytope
    * @return constrained zonotope
    * @ingroup ZonoOpt_SetupFunctions
    *
    * Vpoly is a matrix where each row is a vertex of the polytope.
    */
    std::unique_ptr<ConZono> vrep_2_conzono(const Eigen::Matrix<zono_float, -1, -1>& Vpoly);
} // namespace ZonoOpt


#endif
