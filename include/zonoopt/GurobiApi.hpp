#ifndef ZONOOPT_GUROBI_API_HPP_
#define ZONOOPT_GUROBI_API_HPP_

#ifdef _WIN32
    #include <windows.h>
    typedef HMODULE ZonoOptLibHandle;
    #define ZONOOPT_GET_SYMBOL GetProcAddress
#else
    #include <dlfcn.h>
    typedef void* ZonoOptLibHandle;
    #define ZONOOPT_GET_SYMBOL dlsym
#endif

#include <vector>
#include <string>
#include <memory>

namespace ZonoOpt
{
namespace detail
{

class GurobiApi {
public:

    struct ModelDeleter {
        void operator()(void* model) const;
    };

    using Model = std::unique_ptr<void, ModelDeleter>;
    using Env = std::shared_ptr<void>;

    using LibPtr = std::shared_ptr<void>;

    // function signatures
    typedef int (*GRBloadenv_t)(void**, const char*);
    typedef int (*GRBnewmodel_t)(void*, void**, const char*, int, double*, double*, double*, char*, const char**);
    typedef void (*GRBfreeenv_t)(void*);
    typedef void (*GRBfreemodel_t)(void*);
    typedef int (*GRBaddconstr_t)(void*, int, int*, double*, char, double, const char*);
    typedef int (*GRBoptimize_t)(void*);
    typedef int (*GRBaddqpterms_t)(void*, int, int*, int*, double*);
    typedef int (*GRBupdatemodel_t)(void*);
    typedef int (*GRBgetintattr_t)(void*, const char*, int*);
    typedef int (*GRBgetdblattr_t)(void*, const char*, double*);
    typedef int (*GRBgetdblattrarray_t)(void*, const char*, int, int, double*);
    typedef int (*GRBsetintparam_t)(void*, const char*, int);
    typedef int (*GRBsetdblparam_t)(void*, const char*, double);
    typedef int (*GRBsetstrparam_t)(void*, const char*, const char*);
    typedef void* (*GRBgetenv_t)(void*);
    typedef int (*GRBemptyenv_t)(void**);
    typedef int (*GRBstartenv_t)(void*);

    // function pointers
    GRBloadenv_t GRBloadenv = nullptr;
    GRBnewmodel_t GRBnewmodel = nullptr;
    GRBfreeenv_t GRBfreeenv = nullptr;
    GRBfreemodel_t GRBfreemodel = nullptr;
    GRBaddconstr_t GRBaddconstr = nullptr;
    GRBoptimize_t GRBoptimize = nullptr;
    GRBaddqpterms_t GRBaddqpterms = nullptr;
    GRBupdatemodel_t GRBupdatemodel = nullptr;
    GRBgetintattr_t GRBgetintattr = nullptr;
    GRBgetdblattr_t GRBgetdblattr = nullptr;
    GRBgetdblattrarray_t GRBgetdblattrarray = nullptr;
    GRBsetintparam_t GRBsetintparam = nullptr;
    GRBsetdblparam_t GRBsetdblparam = nullptr;
    GRBsetstrparam_t GRBsetstrparam = nullptr;
    GRBgetenv_t GRBgetenv = nullptr;
    GRBemptyenv_t GRBemptyenv = nullptr;
    GRBstartenv_t GRBstartenv = nullptr;

    // API management
    static GurobiApi& instance();

    bool is_available() const { return _lib_ptr != nullptr; }

    // model and environment
    Env create_env(const std::string& log = "");
    Model create_model(Env& env, const std::string& name, int numvars = 0, double* obj = nullptr,
                       double* lb = nullptr, double* ub = nullptr, char* vtype = nullptr,
                       const char** varnames = nullptr);

    // function wrappers (throw std::runtime_error on failure)
    void add_constr(Model& model, int nnz, int* ind, double* val, char sense, double rhs, const char* name = nullptr);
    void optimize(Model& model);
    void add_qp_terms(Model& model, int nnz, int* qrow, int* qcol, double* qval);
    void update_model(Model& model);
    void get_int_attr(Model& model, const std::string& attr_name, int& value_out);
    void get_dbl_attr(Model& model, const std::string& attr_name, double& value_out);
    void get_dbl_attr_array(Model& model, const std::string& attr_name, int start, int len, double* arr_out);
    void set_int_param(Model& model, const std::string& param_name, int value);
    void set_dbl_param(Model& model, const std::string& param_name, double value);
    void set_str_param(Model& model, const std::string& param_name, const std::string& value);

private:
    LibPtr _lib_ptr = nullptr;
    GurobiApi();
    ~GurobiApi() = default;
    GurobiApi(const GurobiApi&) = delete;
    GurobiApi& operator=(const GurobiApi&) = delete;

    ZonoOptLibHandle load_library();
    ZonoOptLibHandle load_library_from_path(const std::string& path, const std::vector<std::string>& names);
};

} // namespace detail
} // namespace ZonoOpt

#endif
