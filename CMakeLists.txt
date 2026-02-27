cmake_minimum_required(VERSION 3.15...3.27)

option(ZONOOPT_BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)
message(STATUS "ZONOOPT_BUILD_PYTHON_BINDINGS set to: ${ZONOOPT_BUILD_PYTHON_BINDINGS}")

option(ZONOOPT_INSTALL "Build with install rules" OFF)
if(ZONOOPT_BUILD_PYTHON_BINDINGS)
    set(ZONOOPT_INSTALL OFF CACHE BOOL "Disable install build when building Python bindings" FORCE)
endif()
message(STATUS "ZONOOPT_INSTALL set to: ${ZONOOPT_INSTALL}")

option(ZONOOPT_USE_NATIVE "Enable native vectorization (march=native)" OFF)
if(ZONOOPT_INSTALL OR "$ENV{CIBUILDWHEEL}" STREQUAL "1")
    set(ZONOOPT_USE_NATIVE OFF CACHE BOOL "Disable march=native for install and PyPI builds" FORCE)
endif()
message(STATUS "ZONOOPT_USE_NATIVE set to: ${ZONOOPT_USE_NATIVE}")

set(ZONO_FLOAT "double" CACHE STRING "The underlying floating-point type for ZonoOpt library.")
message(STATUS "ZONO_FLOAT set to: ${ZONO_FLOAT}")

set(ZONO_EPS "1e-12" CACHE STRING "Numerical precision for floating point comparisons in ZonoOpt.")
message(STATUS "ZONO_EPS set to: ${ZONO_EPS}")



# project name
if (ZONOOPT_BUILD_PYTHON_BINDINGS)
    # Scikit-build-core sets these values for you, or you can just hard-code the
    # name and version.
    project(
            ${SKBUILD_PROJECT_NAME}
            VERSION ${SKBUILD_PROJECT_VERSION}
            LANGUAGES CXX
    )

    set(TARGET_VAR _core)
else()
    project(ZonoOpt
            VERSION 2.0.1
            LANGUAGES CXX
    )
    set(TARGET_VAR ZonoOpt)
endif()

include(CTest)

# Eigen dependency
include(FetchContent)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if (NOT Eigen3_FOUND)
    if (ZONOOPT_INSTALL)
        message(FATAL_ERROR "Eigen3 not found on system. Required for ZONOOPT_INSTALL.")
    else()
        message(STATUS "Eigen3 not found, fetching 5.0.0 via FetchContent...")
        FetchContent_Declare(
                eigen_default
                GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
                GIT_TAG        "5.0.0"
        )
        FetchContent_MakeAvailable(eigen_default)
    endif()
endif()

if (NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "Failed to fetch Eigen dependency.")
endif()


function(fetch_boost_lib name repo)
  FetchContent_Declare(
    boost_${name}
    GIT_REPOSITORY https://github.com/boostorg/${repo}.git
    GIT_TAG boost-1.90.0
  )
  FetchContent_MakeAvailable(boost_${name})
  list(APPEND BOOST_INCLUDE_DIRS ${boost_${name}_SOURCE_DIR}/include)
  set(BOOST_INCLUDE_DIRS ${BOOST_INCLUDE_DIRS} PARENT_SCOPE)
endfunction()

find_package(Boost CONFIG)
if (NOT Boost_FOUND)
   if (ZONOOPT_INSTALL)
        message(FATAL_ERROR "Boost not found on system. Required for ZONOOPT_INSTALL.")
   else()
        message(STATUS "Boost not found, fetching 1.90.0 via FetchContent...")

        # Fetch interval and its dependencies
        fetch_boost_lib(interval interval)
        fetch_boost_lib(config config)
        fetch_boost_lib(detail detail)
        fetch_boost_lib(assert assert)
        fetch_boost_lib(core core)
        fetch_boost_lib(static_assert static_assert)
        fetch_boost_lib(throw_exception throw_exception)
        fetch_boost_lib(type_traits type_traits)
        fetch_boost_lib(preprocessor preprocessor)
        fetch_boost_lib(concept_check concept_check)
        fetch_boost_lib(numeric_conversion numeric_conversion)

        add_library(boost_interval INTERFACE)
        target_include_directories(boost_interval INTERFACE ${BOOST_INCLUDE_DIRS})
        set(boost_handle boost_interval)
   endif()
else()
    set(boost_handle Boost::boost)
endif()

if (NOT TARGET ${boost_handle})
    message(FATAL_ERROR "Failed to fetch Boost dependency.")
endif()

# python dependencies
if (ZONOOPT_BUILD_PYTHON_BINDINGS)
    # Find the module development requirements (requires FindPython from 3.17 or
    # scikit-build-core's built-in backport)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)
endif()

# source files
set(SOURCE_FILES
    src/SetOperations.cpp
    src/ADMM.cpp
    src/CholeskyUtilities.cpp
    src/HybZono.cpp
    src/ConZono.cpp
    src/Zono.cpp
    src/Point.cpp
    src/EmptySet.cpp
    src/Intervals.cpp
    src/BranchAndBound.cpp
)

# build
if (ZONOOPT_BUILD_PYTHON_BINDINGS)
    # Add a library using FindPython's tooling (pybind11 also provides a helper like
    # this)
    python_add_library(${TARGET_VAR} MODULE
            python/src/zonoopt_py.cpp
            ${SOURCE_FILES}
            WITH_SOABI
    )

    # link libraries
    target_link_libraries(${TARGET_VAR} PRIVATE pybind11::headers)

    # This is passing in the version as a define just as an example
    target_compile_definitions(${TARGET_VAR} PRIVATE VERSION_INFO=${PROJECT_VERSION})

    # The install directory is the output (wheel) directory
    install(TARGETS ${TARGET_VAR} DESTINATION zonoopt)
else()
    add_library(${TARGET_VAR} STATIC
            ${SOURCE_FILES}
    )
endif()

# link Eigen and Boost
target_link_libraries(${TARGET_VAR} PUBLIC Eigen3::Eigen)
target_link_libraries(${TARGET_VAR} PUBLIC ${boost_handle})

# include directories
target_include_directories(${TARGET_VAR} PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# C++ standard
target_compile_features(${TARGET_VAR}
        PUBLIC cxx_std_17
)

# preprocessor definitions
target_compile_definitions(${TARGET_VAR} PRIVATE
        zono_float=${ZONO_FLOAT}
        zono_eps=${ZONO_EPS}
)

# debug flags
if (MSVC)
    target_compile_options(${TARGET_VAR} PRIVATE /bigobj /W3)
else()
    target_compile_options(${TARGET_VAR} PRIVATE -Wall -Wextra -Wpedantic)
endif()

# march native
if (ZONOOPT_USE_NATIVE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TARGET_VAR} PUBLIC -march=native)
    elseif(MSVC)
        target_compile_options(${TARGET_VAR} PUBLIC /arch:AVX2)
    endif()
endif()

# testing
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# install
if (ZONOOPT_INSTALL)
    include(CMakePackageConfigHelpers)

    # Install the target
    install(TARGETS ${TARGET_VAR}
            EXPORT ZonoOptTargets
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include
    )

    # Install headers
    install(DIRECTORY include/ DESTINATION include)

    # Export the targets
    install(EXPORT ZonoOptTargets
            FILE ZonoOptTargets.cmake
            DESTINATION lib/cmake/ZonoOpt
    )

    # Create and install the config files
    configure_package_config_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ZonoOptConfig.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/ZonoOptConfig.cmake"
            INSTALL_DESTINATION lib/cmake/ZonoOpt
    )

    write_basic_package_version_file(
            "${CMAKE_CURRENT_BINARY_DIR}/ZonoOptConfigVersion.cmake"
            VERSION ${PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
    )

    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/ZonoOptConfig.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/ZonoOptConfigVersion.cmake"
            DESTINATION lib/cmake/ZonoOpt
    )
endif()