#include "zonoopt/SCIPApi.hpp"

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace ZonoOpt
{
namespace detail
{

// ---------- Scip RAII wrapper -----------------------------------------------

SCIPApi::Scip::Scip(ScipPtr p, LibPtr lib) : _ptr(p), _lib(std::move(lib)) {}

SCIPApi::Scip::~Scip()
{
    if (_ptr)
    {
        // Use the live function pointer captured from the API singleton.
        auto& api = SCIPApi::instance();
        if (api.is_available() && api.SCIPfree)
        {
            api.SCIPfree(&_ptr);
        }
    }
}

SCIPApi::Scip::Scip(Scip&& other) noexcept : _ptr(other._ptr), _lib(std::move(other._lib))
{
    other._ptr = nullptr;
}

SCIPApi::Scip& SCIPApi::Scip::operator=(Scip&& other) noexcept
{
    if (this != &other)
    {
        if (_ptr)
        {
            auto& api = SCIPApi::instance();
            if (api.is_available() && api.SCIPfree) api.SCIPfree(&_ptr);
        }
        _ptr = other._ptr;
        _lib = std::move(other._lib);
        other._ptr = nullptr;
    }
    return *this;
}

// ---------- SCIPApi singleton ----------------------------------------------

SCIPApi& SCIPApi::instance()
{
    static SCIPApi inst;
    return inst;
}

SCIPApi::Scip SCIPApi::create_scip()
{
    if (!is_available())
    {
        throw std::runtime_error("SCIP API is not available.");
    }

    ScipPtr raw = nullptr;
    int rc = SCIPcreate(&raw);
    if (rc != SCIP_OKAY || !raw)
    {
        throw std::runtime_error("SCIPcreate failed (rc=" + std::to_string(rc) + ").");
    }

    Scip scip(raw, _lib_ptr);

    rc = SCIPincludeDefaultPlugins(raw);
    if (rc != SCIP_OKAY)
    {
        throw std::runtime_error("SCIPincludeDefaultPlugins failed (rc=" + std::to_string(rc) + ").");
    }
    return scip;
}

SCIPApi::SCIPApi()
{
    ZonoOptScipLibHandle h = load_library();
    if (!h)
    {
        _unavailable_reason = "SCIP shared library could not be located (try setting SCIPOPTDIR).";
        return;
    }

    auto deleter = [](void* p) {
        if (!p) return;
#ifdef _WIN32
        FreeLibrary(static_cast<ZonoOptScipLibHandle>(p));
#else
        dlclose(p);
#endif
    };
    _lib_ptr = LibPtr(static_cast<void*>(h), deleter);

    auto sym = [h](const char* name) { return ZONOOPT_SCIP_GET_SYMBOL(h, name); };

    SCIPcreate                          = (SCIPcreate_t)                          sym("SCIPcreate");
    SCIPfree                            = (SCIPfree_t)                            sym("SCIPfree");
    SCIPincludeDefaultPlugins           = (SCIPincludeDefaultPlugins_t)           sym("SCIPincludeDefaultPlugins");
    SCIPcreateProbBasic                 = (SCIPcreateProbBasic_t)                 sym("SCIPcreateProbBasic");
    SCIPfreeProb                        = (SCIPfreeProb_t)                        sym("SCIPfreeProb");
    SCIPsetMessagehdlrQuiet             = (SCIPsetMessagehdlrQuiet_t)             sym("SCIPsetMessagehdlrQuiet");
    SCIPinfinity                        = (SCIPinfinity_t)                        sym("SCIPinfinity");

    SCIPcreateVarBasic                  = (SCIPcreateVarBasic_t)                  sym("SCIPcreateVarBasic");
    SCIPaddVar                          = (SCIPaddVar_t)                          sym("SCIPaddVar");
    SCIPreleaseVar                      = (SCIPreleaseVar_t)                      sym("SCIPreleaseVar");
    SCIPchgVarType                      = (SCIPchgVarType_t)                      sym("SCIPchgVarType");

    SCIPcreateConsBasicLinear           = (SCIPcreateConsBasicLinear_t)           sym("SCIPcreateConsBasicLinear");
    SCIPcreateConsBasicQuadraticNonlinear = (SCIPcreateConsBasicQuadraticNonlinear_t) sym("SCIPcreateConsBasicQuadraticNonlinear");
    SCIPaddCons                         = (SCIPaddCons_t)                         sym("SCIPaddCons");
    SCIPreleaseCons                     = (SCIPreleaseCons_t)                     sym("SCIPreleaseCons");

    SCIPsolve                           = (SCIPsolve_t)                           sym("SCIPsolve");
    SCIPfreeTransform                   = (SCIPfreeTransform_t)                   sym("SCIPfreeTransform");
    SCIPgetStatus                       = (SCIPgetStatus_t)                       sym("SCIPgetStatus");
    SCIPgetBestSol                      = (SCIPgetBestSol_t)                      sym("SCIPgetBestSol");
    SCIPgetNSols                        = (SCIPgetNSols_t)                        sym("SCIPgetNSols");
    SCIPgetSols                         = (SCIPgetSols_t)                         sym("SCIPgetSols");
    SCIPgetSolVal                       = (SCIPgetSolVal_t)                       sym("SCIPgetSolVal");
    SCIPgetSolOrigObj                   = (SCIPgetSolOrigObj_t)                   sym("SCIPgetSolOrigObj");
    SCIPgetSolvingTime                  = (SCIPgetSolvingTime_t)                  sym("SCIPgetSolvingTime");

    SCIPsetBoolParam                    = (SCIPsetBoolParam_t)                    sym("SCIPsetBoolParam");
    SCIPsetIntParam                     = (SCIPsetIntParam_t)                     sym("SCIPsetIntParam");
    SCIPsetLongintParam                 = (SCIPsetLongintParam_t)                 sym("SCIPsetLongintParam");
    SCIPsetRealParam                    = (SCIPsetRealParam_t)                    sym("SCIPsetRealParam");
    SCIPsetCharParam                    = (SCIPsetCharParam_t)                    sym("SCIPsetCharParam");
    SCIPsetStringParam                  = (SCIPsetStringParam_t)                  sym("SCIPsetStringParam");

    SCIPmajorVersion                    = (SCIPmajorVersion_t)                    sym("SCIPmajorVersion");
    SCIPminorVersion                    = (SCIPminorVersion_t)                    sym("SCIPminorVersion");
    SCIPtechVersion                     = (SCIPtechVersion_t)                     sym("SCIPtechVersion");

    // Require all essentials to be present; if any is missing, mark the library as unavailable.
    if (!SCIPcreate || !SCIPfree || !SCIPincludeDefaultPlugins || !SCIPcreateProbBasic ||
        !SCIPcreateVarBasic || !SCIPaddVar || !SCIPreleaseVar ||
        !SCIPcreateConsBasicLinear || !SCIPaddCons || !SCIPreleaseCons ||
        !SCIPsolve || !SCIPgetStatus || !SCIPgetBestSol || !SCIPgetSolVal ||
        !SCIPinfinity || !SCIPmajorVersion)
    {
        _lib_ptr.reset();
        _unavailable_reason = "SCIP shared library is missing one or more required symbols "
                              "(library may be too old or incomplete).";
        return;
    }

    // Determine the loaded SCIP version and adjust enum mappings.
    scip_major = SCIPmajorVersion();
    scip_minor = SCIPminorVersion ? SCIPminorVersion() : 0;
    scip_patch = SCIPtechVersion  ? SCIPtechVersion()  : 0;

    if (scip_major < MIN_SCIP_MAJOR)
    {
        _lib_ptr.reset();
        _unavailable_reason = "Loaded SCIP " + std::to_string(scip_major) + "." +
                              std::to_string(scip_minor) + " is too old; SCIP " +
                              std::to_string(MIN_SCIP_MAJOR) + " or newer is required.";
        return;
    }

    // SCIP_Status enum was reorganized in SCIP 10 (terminating outcomes first).
    // SCIP_STATUS_* constants in SCIPApi.hpp encode the SCIP 10 layout. Older versions
    // were rejected above by the MIN_SCIP_MAJOR check.
}

ZonoOptScipLibHandle SCIPApi::load_library()
{
    // Names to try, ordered from most-specific to most-generic.
    std::vector<std::string> names = {
        "libscip.so", "libscip.so.9", "libscip.so.8", "libscip.so.7",
        "libscipsolver.so",
        "libscip.dylib", "libscip.9.dylib",
        "scip.dll", "libscip.dll"
    };

    // Try SCIPOPTDIR (conventional env var for the SCIP Optimization Suite).
    if (const char* dir = std::getenv("SCIPOPTDIR"))
    {
        const std::string path = std::string(dir) + "/lib/";
        if (auto h = load_library_from_path(path, names)) return h;
    }
    // Try SCIP_DIR as a secondary convention.
    if (const char* dir = std::getenv("SCIP_DIR"))
    {
        const std::string path = std::string(dir) + "/lib/";
        if (auto h = load_library_from_path(path, names)) return h;
    }
    // Fall back to the dynamic linker's default search path.
    return load_library_from_path("", names);
}

ZonoOptScipLibHandle SCIPApi::load_library_from_path(const std::string& path, const std::vector<std::string>& names)
{
    for (const auto& name : names)
    {
        const std::string full = path + name;
        ZonoOptScipLibHandle h = nullptr;
#ifdef _WIN32
        h = LoadLibrary(full.c_str());
#else
        h = dlopen(full.c_str(), RTLD_LAZY);
#endif
        if (h) return h;
    }
    return nullptr;
}

} // namespace detail
} // namespace ZonoOpt
