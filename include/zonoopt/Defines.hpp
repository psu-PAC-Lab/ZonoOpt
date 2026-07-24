#pragma once

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
    constexpr auto pi = static_cast<zono_float>(3.14159265358979323846);
    constexpr auto zero = static_cast<zono_float>(0.0);
    constexpr auto p5 = static_cast<zono_float>(0.5);
    constexpr auto one = static_cast<zono_float>(1.0);
    constexpr auto two = static_cast<zono_float>(2.0);
}