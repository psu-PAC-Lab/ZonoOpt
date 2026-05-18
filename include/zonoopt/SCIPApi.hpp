#ifndef ZONOOPT_SCIP_API_HPP_
#define ZONOOPT_SCIP_API_HPP_

/**
 * @file SCIPApi.hpp
 * @brief Low-level dynamic loader for the SCIP optimization solver.
 *
 * SCIP is loaded at runtime via dlopen / LoadLibrary — there is no build-time
 * dependency on SCIP. If the SCIP shared library cannot be loaded, is_available()
 * returns false and the caller should fall back to the internal solver.
 */

#ifdef _WIN32
    #include <windows.h>
    typedef HMODULE ZonoOptScipLibHandle;
    #define ZONOOPT_SCIP_GET_SYMBOL GetProcAddress
#else
    #include <dlfcn.h>
    typedef void* ZonoOptScipLibHandle;
    #define ZONOOPT_SCIP_GET_SYMBOL dlsym
#endif

#include <memory>
#include <string>
#include <vector>

namespace ZonoOpt
{
namespace detail
{

class SCIPApi {
public:

    // SCIP_RETCODE values that we use (mirrored from scip_retcode.h so we don't need SCIP headers).
    static constexpr int SCIP_OKAY = 1;

    // SCIP_VARTYPE values
    static constexpr int SCIP_VARTYPE_BINARY     = 0;
    static constexpr int SCIP_VARTYPE_INTEGER    = 1;
    static constexpr int SCIP_VARTYPE_IMPLINT    = 2;
    static constexpr int SCIP_VARTYPE_CONTINUOUS = 3;

    // SCIP_Status enum values. Only SCIP 10+ is supported; SCIP 10 reorganized the enum
    // so terminating outcomes come first (OPTIMAL=1, INFEASIBLE=2, UNBOUNDED=3, INFORUNBD=4).
    // SCIP 8 and 9 used different layouts which we have not verified, so they are rejected
    // at load time.
    static constexpr int SCIP_STATUS_OPTIMAL    = 1;
    static constexpr int SCIP_STATUS_INFEASIBLE = 2;
    static constexpr int SCIP_STATUS_UNBOUNDED  = 3;
    static constexpr int SCIP_STATUS_INFORUNBD  = 4;

    // Loaded SCIP version. Zero if no library was successfully loaded.
    int scip_major = 0;
    int scip_minor = 0;
    int scip_patch = 0;

    /// Minimum supported SCIP major version.
    static constexpr int MIN_SCIP_MAJOR = 10;

    // Opaque SCIP_* pointer typedefs (we treat them as void*).
    using ScipPtr  = void*;       // SCIP*
    using VarPtr   = void*;       // SCIP_VAR*
    using ConsPtr  = void*;       // SCIP_CONS*
    using SolPtr   = void*;       // SCIP_SOL*

    using LibPtr = std::shared_ptr<void>;

    // ---- function pointer typedefs ------------------------------------------
    typedef int    (*SCIPcreate_t)(ScipPtr*);
    typedef int    (*SCIPfree_t)(ScipPtr*);
    typedef int    (*SCIPincludeDefaultPlugins_t)(ScipPtr);
    typedef int    (*SCIPcreateProbBasic_t)(ScipPtr, const char*);
    typedef int    (*SCIPfreeProb_t)(ScipPtr);
    typedef int    (*SCIPsetMessagehdlrQuiet_t)(ScipPtr, unsigned int);
    typedef double (*SCIPinfinity_t)(ScipPtr);

    typedef int    (*SCIPcreateVarBasic_t)(ScipPtr, VarPtr*, const char*, double, double, double, int);
    typedef int    (*SCIPaddVar_t)(ScipPtr, VarPtr);
    typedef int    (*SCIPreleaseVar_t)(ScipPtr, VarPtr*);
    typedef int    (*SCIPchgVarType_t)(ScipPtr, VarPtr, int, unsigned int*);

    typedef int    (*SCIPcreateConsBasicLinear_t)(ScipPtr, ConsPtr*, const char*, int, VarPtr*, double*, double, double);
    typedef int    (*SCIPcreateConsBasicQuadraticNonlinear_t)(ScipPtr, ConsPtr*, const char*,
                                                              int, VarPtr*, double*,
                                                              int, VarPtr*, VarPtr*, double*,
                                                              double, double);
    typedef int    (*SCIPaddCons_t)(ScipPtr, ConsPtr);
    typedef int    (*SCIPreleaseCons_t)(ScipPtr, ConsPtr*);

    typedef int    (*SCIPsolve_t)(ScipPtr);
    typedef int    (*SCIPfreeTransform_t)(ScipPtr);
    typedef long long (*SCIPgetNTotalNodes_t)(ScipPtr);
    typedef double (*SCIPgetGap_t)(ScipPtr);
    typedef double (*SCIPgetDualbound_t)(ScipPtr);
    typedef int    (*SCIPgetStatus_t)(ScipPtr);
    typedef SolPtr (*SCIPgetBestSol_t)(ScipPtr);
    typedef int    (*SCIPgetNSols_t)(ScipPtr);
    typedef SolPtr*(*SCIPgetSols_t)(ScipPtr);
    typedef double (*SCIPgetSolVal_t)(ScipPtr, SolPtr, VarPtr);
    typedef double (*SCIPgetSolOrigObj_t)(ScipPtr, SolPtr);
    typedef double (*SCIPgetSolvingTime_t)(ScipPtr);

    typedef int    (*SCIPsetBoolParam_t)(ScipPtr, const char*, unsigned int);
    typedef int    (*SCIPsetIntParam_t)(ScipPtr, const char*, int);
    typedef int    (*SCIPsetLongintParam_t)(ScipPtr, const char*, long long);
    typedef int    (*SCIPsetRealParam_t)(ScipPtr, const char*, double);
    typedef int    (*SCIPsetCharParam_t)(ScipPtr, const char*, char);
    typedef int    (*SCIPsetStringParam_t)(ScipPtr, const char*, const char*);

    typedef int    (*SCIPmajorVersion_t)();
    typedef int    (*SCIPminorVersion_t)();
    typedef int    (*SCIPtechVersion_t)();

    // ---- function pointers --------------------------------------------------
    SCIPcreate_t                          SCIPcreate                          = nullptr;
    SCIPfree_t                            SCIPfree                            = nullptr;
    SCIPincludeDefaultPlugins_t           SCIPincludeDefaultPlugins           = nullptr;
    SCIPcreateProbBasic_t                 SCIPcreateProbBasic                 = nullptr;
    SCIPfreeProb_t                        SCIPfreeProb                        = nullptr;
    SCIPsetMessagehdlrQuiet_t             SCIPsetMessagehdlrQuiet             = nullptr;
    SCIPinfinity_t                        SCIPinfinity                        = nullptr;

    SCIPcreateVarBasic_t                  SCIPcreateVarBasic                  = nullptr;
    SCIPaddVar_t                          SCIPaddVar                          = nullptr;
    SCIPreleaseVar_t                      SCIPreleaseVar                      = nullptr;
    SCIPchgVarType_t                      SCIPchgVarType                      = nullptr;

    SCIPcreateConsBasicLinear_t           SCIPcreateConsBasicLinear           = nullptr;
    SCIPcreateConsBasicQuadraticNonlinear_t SCIPcreateConsBasicQuadraticNonlinear = nullptr;
    SCIPaddCons_t                         SCIPaddCons                         = nullptr;
    SCIPreleaseCons_t                     SCIPreleaseCons                     = nullptr;

    SCIPsolve_t                           SCIPsolve                           = nullptr;
    SCIPfreeTransform_t                   SCIPfreeTransform                   = nullptr;
    SCIPgetNTotalNodes_t                  SCIPgetNTotalNodes                  = nullptr;
    SCIPgetGap_t                          SCIPgetGap                          = nullptr;
    SCIPgetDualbound_t                    SCIPgetDualbound                    = nullptr;
    SCIPgetStatus_t                       SCIPgetStatus                       = nullptr;
    SCIPgetBestSol_t                      SCIPgetBestSol                      = nullptr;
    SCIPgetNSols_t                        SCIPgetNSols                        = nullptr;
    SCIPgetSols_t                         SCIPgetSols                         = nullptr;
    SCIPgetSolVal_t                       SCIPgetSolVal                       = nullptr;
    SCIPgetSolOrigObj_t                   SCIPgetSolOrigObj                   = nullptr;
    SCIPgetSolvingTime_t                  SCIPgetSolvingTime                  = nullptr;

    SCIPsetBoolParam_t                    SCIPsetBoolParam                    = nullptr;
    SCIPsetIntParam_t                     SCIPsetIntParam                     = nullptr;
    SCIPsetLongintParam_t                 SCIPsetLongintParam                 = nullptr;
    SCIPsetRealParam_t                    SCIPsetRealParam                    = nullptr;
    SCIPsetCharParam_t                    SCIPsetCharParam                    = nullptr;
    SCIPsetStringParam_t                  SCIPsetStringParam                  = nullptr;

    SCIPmajorVersion_t                    SCIPmajorVersion                    = nullptr;
    SCIPminorVersion_t                    SCIPminorVersion                    = nullptr;
    SCIPtechVersion_t                     SCIPtechVersion                     = nullptr;

    // ---- API management -----------------------------------------------------
    static SCIPApi& instance();

    /// True if SCIP was successfully loaded and all essential symbols resolved.
    bool is_available() const { return _lib_ptr != nullptr; }

    /// If is_available() is false, returns a short explanation of why
    /// (e.g., "library not found", "SCIP version too old"). Empty when available.
    const std::string& unavailable_reason() const { return _unavailable_reason; }

    /**
     * @brief RAII wrapper around a SCIP* environment.
     *
     * Owns the SCIP*; destructor calls SCIPfree. Also keeps a copy of the library
     * handle so the shared library stays loaded until the last Scip is destroyed.
     */
    class Scip
    {
    public:
        Scip() = default;
        Scip(ScipPtr p, LibPtr lib);
        ~Scip();
        Scip(const Scip&) = delete;
        Scip& operator=(const Scip&) = delete;
        Scip(Scip&& other) noexcept;
        Scip& operator=(Scip&& other) noexcept;

        ScipPtr get() const { return _ptr; }
        explicit operator bool() const { return _ptr != nullptr; }

    private:
        ScipPtr _ptr = nullptr;
        LibPtr  _lib;
    };

    /// Create a new SCIP environment with default plugins registered.
    Scip create_scip();

private:
    LibPtr _lib_ptr;
    std::string _unavailable_reason;
    SCIPApi();
    ~SCIPApi() = default;
    SCIPApi(const SCIPApi&) = delete;
    SCIPApi& operator=(const SCIPApi&) = delete;

    ZonoOptScipLibHandle load_library();
    ZonoOptScipLibHandle load_library_from_path(const std::string& path, const std::vector<std::string>& names);
};

} // namespace detail
} // namespace ZonoOpt

#endif
