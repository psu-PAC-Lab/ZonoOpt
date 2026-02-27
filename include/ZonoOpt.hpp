#ifndef ZONOOPT_HPP_
#define ZONOOPT_HPP_

#include <sstream>
#include <iostream>

/**
 * @file ZonoOpt.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief ZonoOpt library main header file.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 */

/**
 * @brief Set operations for ZonoOpt library
 * @defgroup ZonoOpt_SetOperations Set Operations
 */

/**
 * @brief Setup functions for ZonoOpt library
 * @defgroup ZonoOpt_SetupFunctions Setup Functions
 */

/**
 * @brief Preprocessor macros for ZonoOpt library 
 * @defgroup ZonoOpt_PreprocessorMacros Preprocessor Macros
 */

/**
 * @brief Typedefs for ZonoOpt library
 * @defgroup ZonoOpt_Typedefs Typedefs
 */

#define EIGEN_MPL2_ONLY // Disable features licensed under LGPL

// ZonoOpt preprocessor directives
/**
 * @brief Defines the floating-point type used in ZonoOpt.
 * @ingroup ZonoOpt_PreprocessorMacros
 */
#ifndef zono_float
    #define zono_float double
#endif

/**
 * @brief Defines the precision used for floating point comparisons in ZonoOpt.
 * @ingroup ZonoOpt_PreprocessorMacros
 */
#ifndef zono_eps
    #define zono_eps Eigen::NumTraits<zono_float>::dummy_precision()
#endif

// constants
namespace ZonoOpt::detail
{
    constexpr zono_float pi = static_cast<zono_float>(3.14159265358979323846);
    constexpr zono_float zero = static_cast<zono_float>(0.0);
    constexpr zono_float p5 = static_cast<zono_float>(0.5);
    constexpr zono_float one = static_cast<zono_float>(1.0);
    constexpr zono_float two = static_cast<zono_float>(2.0);
}

// includes
#include "zonoopt/ADMM.hpp"
#include "zonoopt/CholeskyUtilities.hpp"
#include "zonoopt/ConZono.hpp"
#include "zonoopt/EmptySet.hpp"
#include "zonoopt/GenUtilities.hpp"
#include "zonoopt/HybZono.hpp"
#include "zonoopt/Inequality.hpp"
#include "zonoopt/Intervals.hpp"
#include "zonoopt/BnbDataStructures.hpp"
#include "zonoopt/BranchAndBound.hpp"
#include "zonoopt/Point.hpp"
#include "zonoopt/SolverDataStructures.hpp"
#include "zonoopt/SparseMatrixUtilities.hpp"
#include "zonoopt/Zono.hpp"

// typedef
namespace ZonoOpt
{
    /**
     * @brief Type alias for a unique pointer to a (polymorphic) HybZono object.
     * @ingroup ZonoOpt_Typedefs    
     */
    typedef std::unique_ptr<HybZono> ZonoPtr;
}

#endif