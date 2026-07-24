#pragma once

#include <memory>

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

// includes
#include "zonoopt/Point.hpp"
#include "zonoopt/Zono.hpp"
#include "zonoopt/ConZono.hpp"
#include "zonoopt/HybZono.hpp"
#include "zonoopt/EmptySet.hpp"
#include "zonoopt/Interval.hpp"
#include "zonoopt/Box.hpp"
#include "zonoopt/IntervalMatrix.hpp"
#include "zonoopt/SolverDataStructures.hpp"
#include "zonoopt/GurobiSettings.hpp"
#include "zonoopt/SCIPSettings.hpp"
#include "zonoopt/ZonoJson.hpp"

// typedef
namespace ZonoOpt
{
    /**
     * @brief Type alias for a unique pointer to a (polymorphic) HybZono object.
     * @ingroup ZonoOpt_Typedefs    
     */
    typedef std::unique_ptr<HybZono> ZonoPtr;
}