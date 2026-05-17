#include "zonoopt/GurobiApi.hpp"

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace ZonoOpt
{
namespace detail
{

void GurobiApi::ModelDeleter::operator()(void* model) const
{
    auto& api = GurobiApi::instance();
    if (api.is_available() && model)
    {
        api.GRBfreemodel(model);
    }
}

GurobiApi& GurobiApi::instance()
{
    static GurobiApi inst;
    return inst;
}

GurobiApi::Model GurobiApi::create_model(Env& env, const std::string& name, int numvars,
                                         double* obj, double* lb, double* ub, char* vtype,
                                         const char** varnames)
{
    void* raw_ptr = nullptr;

    int err = this->GRBnewmodel(env.get(), &raw_ptr, name.c_str(), numvars, obj, lb, ub, vtype, varnames);

    if (err || !raw_ptr)
    {
        throw std::runtime_error("Failed to create Gurobi model. Error code: " + std::to_string(err));
    }

    return Model(raw_ptr);
}

GurobiApi::Env GurobiApi::create_env(const std::string& log)
{
    void* raw = nullptr;
    const char* log_arg = log.empty() ? nullptr : log.c_str();
    int err = GRBloadenv(&raw, log_arg);

    if (err || !raw)
    {
        throw std::runtime_error("Failed to create Gurobi environment. Error code: " + std::to_string(err));
    }

    // Capture a copy of the library handle so the .so stays loaded
    // until this specific environment is destroyed.
    auto lib_handle = this->_lib_ptr;
    auto free_func = this->GRBfreeenv;

    return Env(raw, [lib_handle, free_func](void* ptr) {
        if (ptr && free_func)
        {
            free_func(ptr);
        }
    });
}

void GurobiApi::add_constr(Model& model, int nnz, int* ind, double* val, char sense, double rhs, const char* name)
{
    if (!model || !this->GRBaddconstr)
    {
        throw std::runtime_error("Model is null or GRBaddconstr not available.");
    }

    int err = this->GRBaddconstr(model.get(), nnz, ind, val, sense, rhs, name);

    if (err)
    {
        throw std::runtime_error("Failed to add constraint. Error code: " + std::to_string(err));
    }
}

void GurobiApi::optimize(Model& model)
{
    if (!model || !this->GRBoptimize)
    {
        throw std::runtime_error("Model is null or GRBoptimize not available.");
    }

    int err = this->GRBoptimize(model.get());

    if (err)
    {
        throw std::runtime_error("Failed to optimize model. Error code: " + std::to_string(err));
    }
}

void GurobiApi::get_int_attr(Model& model, const std::string& attr_name, int& value_out)
{
    if (!model || !this->GRBgetintattr)
    {
        throw std::runtime_error("Model is null or GRBgetintattr not available.");
    }

    int err = this->GRBgetintattr(model.get(), attr_name.c_str(), &value_out);

    if (err)
    {
        throw std::runtime_error("Failed to get integer attribute. Error code: " + std::to_string(err));
    }
}

void GurobiApi::get_dbl_attr(Model& model, const std::string& attr_name, double& value_out)
{
    if (!model || !this->GRBgetdblattr)
    {
        throw std::runtime_error("Model is null or GRBgetdblattr not available.");
    }

    int err = this->GRBgetdblattr(model.get(), attr_name.c_str(), &value_out);

    if (err)
    {
        throw std::runtime_error("Failed to get double attribute. Error code: " + std::to_string(err));
    }
}

void GurobiApi::get_dbl_attr_array(Model& model, const std::string& attr_name, int start, int len, double* arr_out)
{
    if (!model || !this->GRBgetdblattrarray)
    {
        throw std::runtime_error("Model is null or GRBgetdblattrarray not available.");
    }

    int err = this->GRBgetdblattrarray(model.get(), attr_name.c_str(), start, len, arr_out);

    if (err)
    {
        throw std::runtime_error("Failed to get attribute array. Error code: " + std::to_string(err));
    }
}

void GurobiApi::add_qp_terms(Model& model, int nnz, int* qrow, int* qcol, double* qval)
{
    if (!model || !this->GRBaddqpterms)
    {
        throw std::runtime_error("Model is null or GRBaddqpterms not available.");
    }

    int err = this->GRBaddqpterms(model.get(), nnz, qrow, qcol, qval);

    if (err)
    {
        throw std::runtime_error("Failed to add quadratic terms. Error code: " + std::to_string(err));
    }
}

void GurobiApi::update_model(Model& model)
{
    if (!model || !this->GRBupdatemodel)
    {
        throw std::runtime_error("Model is null or GRBupdatemodel not available.");
    }

    int err = this->GRBupdatemodel(model.get());

    if (err)
    {
        throw std::runtime_error("Failed to update model. Error code: " + std::to_string(err));
    }
}

void GurobiApi::set_int_param(Model& model, const std::string& param_name, int value)
{
    if (!model || !this->GRBsetintparam || !this->GRBgetenv)
    {
        throw std::runtime_error("Model is null or GRBsetintparam/GRBgetenv not available.");
    }

    void* model_env = this->GRBgetenv(model.get());
    if (!model_env)
    {
        throw std::runtime_error("Failed to get model environment.");
    }

    int err = this->GRBsetintparam(model_env, param_name.c_str(), value);

    if (err)
    {
        throw std::runtime_error("Failed to set integer parameter '" + param_name + "'. Error code: " +
                                 std::to_string(err));
    }
}

void GurobiApi::set_dbl_param(Model& model, const std::string& param_name, double value)
{
    if (!model || !this->GRBsetdblparam || !this->GRBgetenv)
    {
        throw std::runtime_error("Model is null or GRBsetdblparam/GRBgetenv not available.");
    }

    void* model_env = this->GRBgetenv(model.get());
    if (!model_env)
    {
        throw std::runtime_error("Failed to get model environment.");
    }

    int err = this->GRBsetdblparam(model_env, param_name.c_str(), value);

    if (err)
    {
        throw std::runtime_error("Failed to set double parameter '" + param_name + "'. Error code: " +
                                 std::to_string(err));
    }
}

void GurobiApi::set_str_param(Model& model, const std::string& param_name, const std::string& value)
{
    if (!model || !this->GRBsetstrparam || !this->GRBgetenv)
    {
        throw std::runtime_error("Model is null or GRBsetstrparam/GRBgetenv not available.");
    }

    void* model_env = this->GRBgetenv(model.get());
    if (!model_env)
    {
        throw std::runtime_error("Failed to get model environment.");
    }

    int err = this->GRBsetstrparam(model_env, param_name.c_str(), value.c_str());

    if (err)
    {
        throw std::runtime_error("Failed to set string parameter '" + param_name + "'. Error code: " +
                                 std::to_string(err));
    }
}

GurobiApi::GurobiApi()
{
    ZonoOptLibHandle h = load_library();
    if (h)
    {
        auto deleter = [](void* p) {
            if (!p) return;
#ifdef _WIN32
            FreeLibrary(static_cast<ZonoOptLibHandle>(p));
#else
            dlclose(p);
#endif
        };

        _lib_ptr = LibPtr(static_cast<void*>(h), deleter);

        GRBloadenv         = (GRBloadenv_t)         ZONOOPT_GET_SYMBOL(h, "GRBloadenv");
        GRBnewmodel        = (GRBnewmodel_t)        ZONOOPT_GET_SYMBOL(h, "GRBnewmodel");
        GRBfreeenv         = (GRBfreeenv_t)         ZONOOPT_GET_SYMBOL(h, "GRBfreeenv");
        GRBfreemodel       = (GRBfreemodel_t)       ZONOOPT_GET_SYMBOL(h, "GRBfreemodel");
        GRBaddconstr       = (GRBaddconstr_t)       ZONOOPT_GET_SYMBOL(h, "GRBaddconstr");
        GRBoptimize        = (GRBoptimize_t)        ZONOOPT_GET_SYMBOL(h, "GRBoptimize");
        GRBgetdblattrarray = (GRBgetdblattrarray_t) ZONOOPT_GET_SYMBOL(h, "GRBgetdblattrarray");
        GRBaddqpterms      = (GRBaddqpterms_t)      ZONOOPT_GET_SYMBOL(h, "GRBaddqpterms");
        GRBupdatemodel     = (GRBupdatemodel_t)     ZONOOPT_GET_SYMBOL(h, "GRBupdatemodel");
        GRBgetintattr      = (GRBgetintattr_t)      ZONOOPT_GET_SYMBOL(h, "GRBgetintattr");
        GRBgetdblattr      = (GRBgetdblattr_t)      ZONOOPT_GET_SYMBOL(h, "GRBgetdblattr");
        GRBsetintparam     = (GRBsetintparam_t)     ZONOOPT_GET_SYMBOL(h, "GRBsetintparam");
        GRBsetdblparam     = (GRBsetdblparam_t)     ZONOOPT_GET_SYMBOL(h, "GRBsetdblparam");
        GRBsetstrparam     = (GRBsetstrparam_t)     ZONOOPT_GET_SYMBOL(h, "GRBsetstrparam");
        GRBgetenv          = (GRBgetenv_t)          ZONOOPT_GET_SYMBOL(h, "GRBgetenv");
        GRBemptyenv        = (GRBemptyenv_t)        ZONOOPT_GET_SYMBOL(h, "GRBemptyenv");
        GRBstartenv        = (GRBstartenv_t)        ZONOOPT_GET_SYMBOL(h, "GRBstartenv");

        // require all essentials to be resolved; otherwise treat as unavailable
        if (!GRBloadenv || !GRBnewmodel || !GRBfreeenv || !GRBfreemodel || !GRBaddconstr ||
            !GRBoptimize || !GRBgetdblattrarray || !GRBgetintattr)
        {
            _lib_ptr.reset();
        }
    }
}

ZonoOptLibHandle GurobiApi::load_library()
{
    std::vector<std::string> names = {
        "libgurobi130.so", "libgurobi120.so", "libgurobi110.so", "libgurobi100.so", "libgurobi.so",
        "gurobi130.dll", "gurobi120.dll", "gurobi110.dll", "gurobi100.dll", "gurobi.dll",
        "libgurobi130.dylib", "libgurobi120.dylib", "libgurobi110.dylib", "libgurobi.dylib"
    };

    const char* home = std::getenv("GUROBI_HOME");

    if (home)
    {
        const std::string path = std::string(home) + "/lib/";
        ZonoOptLibHandle h = load_library_from_path(path, names);
        if (h) return h;
    }

    return load_library_from_path("", names);
}

ZonoOptLibHandle GurobiApi::load_library_from_path(const std::string& path, const std::vector<std::string>& names)
{
    ZonoOptLibHandle h = nullptr;
    for (const auto& name : names)
    {
        const std::string full_path = path + name;
#ifdef _WIN32
        h = LoadLibrary(full_path.c_str());
#else
        h = dlopen(full_path.c_str(), RTLD_LAZY);
#endif
        if (h) return h;
    }
    return nullptr;
}

} // namespace detail
} // namespace ZonoOpt
