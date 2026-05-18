#include "zonoopt/GurobiSolver.hpp"
#include "zonoopt/GurobiApi.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <vector>
#include <string>
#include <stdexcept>
#include <utility>

namespace ZonoOpt
{
namespace detail
{

namespace
{

// Gurobi status codes (from gurobi_c.h) — replicated here so we don't need Gurobi at compile time.
constexpr int GRB_OPTIMAL         = 2;
constexpr int GRB_INFEASIBLE      = 3;
constexpr int GRB_INF_OR_UNBD     = 4;
constexpr int GRB_SUBOPTIMAL      = 13;

struct GurobiResult
{
    int status = 0;
    double runtime = 0.0;
    double objective = 0.0;
    Eigen::VectorXd y;
    bool got_solution = false;

    // Solver-native metadata fetched after optimization for GurobiSolverResults.
    double iter_count = 0.0;
    double node_count = 0.0;
    double mip_gap = 0.0;
    double obj_bound = 0.0;
};

/**
 * Read Gurobi's solver-native metadata attributes (best-effort).
 * "NodeCount" / "ObjBound" / "MIPGap" are MIP-only and silently skipped on LP/QP models.
 */
void fetch_gurobi_metadata(GurobiApi& api, GurobiApi::Model& model, GurobiResult& gr)
{
    auto try_dbl = [&](const char* name, double& out) {
        try { api.get_dbl_attr(model, name, out); }
        catch (const std::exception&) { /* attribute may not apply to this model type */ }
    };
    try_dbl("IterCount", gr.iter_count);
    try_dbl("NodeCount", gr.node_count);
    try_dbl("MIPGap",    gr.mip_gap);
    try_dbl("ObjBound",  gr.obj_bound);
}

template <typename T>
Eigen::VectorXd to_double_vec(const Eigen::Vector<T, -1>& v)
{
    return v.template cast<double>();
}

template <typename T>
Eigen::SparseMatrix<double, Eigen::RowMajor> to_double_rowmajor(const Eigen::SparseMatrix<T>& M)
{
    if (M.nonZeros() == 0)
    {
        return Eigen::SparseMatrix<double, Eigen::RowMajor>(M.rows(), M.cols());
    }
    Eigen::SparseMatrix<double> Md = M.template cast<double>();
    return Eigen::SparseMatrix<double, Eigen::RowMajor>(Md);
}

/**
 * Build a Gurobi model with the given QP/MIQP data and apply common settings.
 * Returns ownership of the env and model so the caller can drive optimization
 * and post-optimization attribute reads.
 */
struct PreparedModel
{
    GurobiApi::Env env;
    GurobiApi::Model model;
    int n = 0;
};

void apply_gurobi_settings(GurobiApi& api, GurobiApi::Model& model, const GurobiSettings& s)
{
    auto apply_int = [&](const char* name, const std::optional<int>& v)         { if (v) api.set_int_param(model, name, *v); };
    auto apply_dbl = [&](const char* name, const std::optional<double>& v)      { if (v) api.set_dbl_param(model, name, *v); };
    auto apply_str = [&](const char* name, const std::optional<std::string>& v) { if (v) api.set_str_param(model, name, *v); };

    // Termination
    apply_dbl("TimeLimit", s.TimeLimit);
    apply_dbl("WorkLimit", s.WorkLimit);
    apply_dbl("MemLimit",  s.MemLimit);
    apply_int("SolutionLimit", s.SolutionLimit);

    // Tolerances
    apply_dbl("MIPGap",         s.MIPGap);
    apply_dbl("MIPGapAbs",      s.MIPGapAbs);
    apply_dbl("FeasibilityTol", s.FeasibilityTol);
    apply_dbl("OptimalityTol",  s.OptimalityTol);
    apply_dbl("IntFeasTol",     s.IntFeasTol);

    // Algorithm / behavior
    apply_int("Method",       s.Method);
    apply_int("Presolve",     s.Presolve);
    apply_int("Cuts",         s.Cuts);
    apply_int("MIPFocus",     s.MIPFocus);
    apply_int("NumericFocus", s.NumericFocus);
    apply_dbl("Heuristics",   s.Heuristics);
    apply_int("Threads",      s.Threads);
    apply_int("Seed",         s.Seed);

    // Solution pool
    apply_int("PoolSolutions",  s.PoolSolutions);
    apply_int("PoolSearchMode", s.PoolSearchMode);
    apply_dbl("PoolGap",        s.PoolGap);
    apply_dbl("PoolGapAbs",     s.PoolGapAbs);

    // Logging
    apply_int("OutputFlag",   s.OutputFlag);
    apply_int("LogToConsole", s.LogToConsole);
    apply_str("LogFile",      s.LogFile);

    // Escape-hatch maps for any parameter not exposed above.
    for (const auto& kv : s.int_params) api.set_int_param(model, kv.first, kv.second);
    for (const auto& kv : s.dbl_params) api.set_dbl_param(model, kv.first, kv.second);
    for (const auto& kv : s.str_params) api.set_str_param(model, kv.first, kv.second);
}

PreparedModel prepare_model(const Eigen::SparseMatrix<double, Eigen::RowMajor>& P_rm,
                            const Eigen::VectorXd& q_d,
                            const Eigen::SparseMatrix<double, Eigen::RowMajor>& A_rm,
                            const Eigen::VectorXd& b_d,
                            const Eigen::VectorXd& lb_d,
                            const Eigen::VectorXd& ub_d,
                            const std::vector<char>& vtype,
                            const GurobiSettings& settings)
{
    GurobiApi& api = GurobiApi::instance();
    if (!api.is_available())
    {
        throw std::runtime_error("Gurobi API is not available.");
    }

    PreparedModel pm;
    pm.n = static_cast<int>(q_d.size());

    // LogFile is special: GRBloadenv takes the log filename directly so the env logs
    // creation messages too. We also re-apply it via set_str_param below for parity.
    pm.env = api.create_env(settings.LogFile.value_or(""));
    pm.model = api.create_model(pm.env, "zonoopt_model", pm.n,
                                const_cast<double*>(q_d.data()),
                                const_cast<double*>(lb_d.data()),
                                const_cast<double*>(ub_d.data()),
                                const_cast<char*>(vtype.data()));

    apply_gurobi_settings(api, pm.model, settings);

    for (int k = 0; k < A_rm.outerSize(); ++k)
    {
        std::vector<int> ind;
        std::vector<double> val;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_rm, k); it; ++it)
        {
            ind.push_back(static_cast<int>(it.col()));
            val.push_back(it.value());
        }
        if (!ind.empty())
        {
            api.add_constr(pm.model, static_cast<int>(ind.size()), ind.data(), val.data(),
                           '=', b_d(k));
        }
    }

    if (P_rm.nonZeros() > 0)
    {
        std::vector<int> qrow, qcol;
        std::vector<double> qval;
        qrow.reserve(P_rm.nonZeros());
        qcol.reserve(P_rm.nonZeros());
        qval.reserve(P_rm.nonZeros());
        for (int k = 0; k < P_rm.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(P_rm, k); it; ++it)
            {
                if (it.col() > it.row()) continue;
                const double v = (it.row() == it.col()) ? 0.5 * it.value() : it.value();
                qrow.push_back(static_cast<int>(it.row()));
                qcol.push_back(static_cast<int>(it.col()));
                qval.push_back(v);
            }
        }
        if (!qrow.empty())
        {
            api.add_qp_terms(pm.model, static_cast<int>(qrow.size()),
                             qrow.data(), qcol.data(), qval.data());
        }
    }

    return pm;
}

GurobiResult optimize_single(PreparedModel& pm)
{
    GurobiApi& api = GurobiApi::instance();
    GurobiResult result;
    result.y = Eigen::VectorXd::Zero(pm.n);

    auto t_start = std::chrono::steady_clock::now();
    api.optimize(pm.model);
    auto t_end = std::chrono::steady_clock::now();
    result.runtime = std::chrono::duration<double>(t_end - t_start).count();

    api.get_int_attr(pm.model, "Status", result.status);

    if (result.status == GRB_OPTIMAL || result.status == GRB_SUBOPTIMAL)
    {
        try
        {
            api.get_dbl_attr_array(pm.model, "X", 0, pm.n, result.y.data());
            api.get_dbl_attr(pm.model, "ObjVal", result.objective);
            result.got_solution = true;
        }
        catch (const std::exception&)
        {
            result.got_solution = false;
        }
    }

    fetch_gurobi_metadata(api, pm.model, result);

    return result;
}

/**
 * Drive Gurobi's solution pool to enumerate up to n_sols feasible solutions.
 *
 * The caller's request for n_sols overrides whatever PoolSolutions value was already
 * applied from settings. PoolSearchMode is set to 2 (find n best) unless the caller's
 * settings already specified a non-default value via prepare_model.
 */
std::vector<GurobiResult> optimize_pool(PreparedModel& pm, int n_sols, int pool_search_mode)
{
    GurobiApi& api = GurobiApi::instance();
    std::vector<GurobiResult> results;

    if (n_sols < 1) n_sols = 1;

    api.set_int_param(pm.model, "PoolSolutions", n_sols);
    api.set_int_param(pm.model, "PoolSearchMode", pool_search_mode);

    auto t_start = std::chrono::steady_clock::now();
    api.optimize(pm.model);
    auto t_end = std::chrono::steady_clock::now();
    const double runtime = std::chrono::duration<double>(t_end - t_start).count();

    int status = 0;
    api.get_int_attr(pm.model, "Status", status);

    if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD)
    {
        GurobiResult gr;
        gr.status = status;
        gr.runtime = runtime;
        gr.y = Eigen::VectorXd::Zero(pm.n);
        results.push_back(gr);
        return results;
    }

    int sol_count = 0;
    try
    {
        api.get_int_attr(pm.model, "SolCount", sol_count);
    }
    catch (const std::exception&)
    {
        sol_count = 0;
    }

    if (sol_count == 0)
    {
        GurobiResult gr;
        gr.status = status;
        gr.runtime = runtime;
        gr.y = Eigen::VectorXd::Zero(pm.n);
        results.push_back(gr);
        return results;
    }

    // Fetch model-level metadata once — these attributes describe the optimize call
    // overall, not any particular pool entry.
    GurobiResult shared_meta;
    shared_meta.status = status;
    shared_meta.runtime = runtime;
    fetch_gurobi_metadata(api, pm.model, shared_meta);

    for (int k = 0; k < sol_count; ++k)
    {
        GurobiResult gr = shared_meta;  // copy status/runtime/metadata
        gr.y = Eigen::VectorXd::Zero(pm.n);

        try
        {
            api.set_int_param(pm.model, "SolutionNumber", k);
            api.get_dbl_attr_array(pm.model, "Xn", 0, pm.n, gr.y.data());
            api.get_dbl_attr(pm.model, "PoolObjVal", gr.objective);
            gr.got_solution = true;
        }
        catch (const std::exception&)
        {
            gr.got_solution = false;
        }
        results.push_back(gr);
    }

    return results;
}

OptSolution build_opt_solution(const GurobiResult& gr, const Eigen::Vector<zono_float, -1>& z,
                               zono_float constant_offset)
{
    OptSolution sol;
    sol.z = z;
    sol.run_time = gr.runtime;
    sol.iter = 0;

    if (gr.status == GRB_INFEASIBLE || gr.status == GRB_INF_OR_UNBD)
    {
        sol.infeasible = true;
        sol.converged = false;
        sol.J = std::numeric_limits<zono_float>::infinity();
    }
    else if (gr.got_solution)
    {
        sol.infeasible = false;
        sol.converged = (gr.status == GRB_OPTIMAL || gr.status == GRB_SUBOPTIMAL);
        sol.J = static_cast<zono_float>(gr.objective) + constant_offset;
    }
    else
    {
        sol.infeasible = false;
        sol.converged = false;
        sol.J = -std::numeric_limits<zono_float>::infinity();
    }

    // Attach Gurobi-native solution metadata so callers can read raw status, node count, etc.
    auto gres = std::make_shared<GurobiSolverResults>();
    gres->status     = gr.status;
    gres->iter_count = gr.iter_count;
    gres->node_count = gr.node_count;
    gres->mip_gap    = gr.mip_gap;
    gres->obj_bound  = gr.obj_bound;
    sol.external_results = gres;

    return sol;
}

/**
 * Prepare the MIQP data: optionally applies the xi = M*y + s substitution when binary
 * variables are in {-1, 1} form so Gurobi (which uses {0, 1} binaries) can solve in y-space.
 * Outputs the row-major P, A, vectors, vtype, and the constant offset / transform info
 * needed to recover xi from y after solving.
 */
struct MiqpPrep
{
    Eigen::SparseMatrix<double, Eigen::RowMajor> P_rm;
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_rm;
    Eigen::VectorXd q_d;
    Eigen::VectorXd b_d;
    Eigen::VectorXd lb_d;
    Eigen::VectorXd ub_d;
    std::vector<char> vtype;
    double const_shift = 0.0;
    bool needs_transform = false;
    int bin_start = 0;
    int bin_count = 0;

    Eigen::VectorXd recover_xi(const Eigen::VectorXd& y) const
    {
        if (!needs_transform) return y;
        Eigen::VectorXd xi = y;
        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            xi(i) = 2.0 * y(i) - 1.0;
        }
        return xi;
    }
};

MiqpPrep prep_miqp(const Eigen::SparseMatrix<zono_float>& P,
                   const Eigen::Vector<zono_float, -1>& q,
                   const Eigen::SparseMatrix<zono_float>& A,
                   const Eigen::Vector<zono_float, -1>& b,
                   const Eigen::Vector<zono_float, -1>& xi_lb,
                   const Eigen::Vector<zono_float, -1>& xi_ub,
                   int bin_start, int bin_count,
                   bool zero_one_form)
{
    const int n = static_cast<int>(q.size());

    Eigen::SparseMatrix<double> P_d = (P.nonZeros() > 0)
                                        ? Eigen::SparseMatrix<double>(P.template cast<double>())
                                        : Eigen::SparseMatrix<double>(n, n);
    Eigen::SparseMatrix<double> A_d = (A.nonZeros() > 0)
                                        ? Eigen::SparseMatrix<double>(A.template cast<double>())
                                        : Eigen::SparseMatrix<double>(static_cast<int>(b.size()), n);
    Eigen::VectorXd q_d  = to_double_vec(q);
    Eigen::VectorXd b_d  = to_double_vec(b);
    Eigen::VectorXd lb_d = to_double_vec(xi_lb);
    Eigen::VectorXd ub_d = to_double_vec(xi_ub);

    MiqpPrep out;
    out.bin_start = bin_start;
    out.bin_count = bin_count;
    out.needs_transform = (bin_count > 0) && !zero_one_form;

    if (out.needs_transform)
    {
        Eigen::VectorXd m_vec = Eigen::VectorXd::Ones(n);
        Eigen::VectorXd s_vec = Eigen::VectorXd::Zero(n);
        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            m_vec(i) = 2.0;
            s_vec(i) = -1.0;
        }

        Eigen::VectorXd P_s = (P_d.nonZeros() > 0) ? Eigen::VectorXd(P_d * s_vec) : Eigen::VectorXd::Zero(n);
        Eigen::VectorXd new_q = m_vec.asDiagonal() * (P_s + q_d);
        out.const_shift = 0.5 * s_vec.dot(P_s) + q_d.dot(s_vec);

        if (P_d.nonZeros() > 0)
        {
            P_d = m_vec.asDiagonal() * P_d * m_vec.asDiagonal();
        }
        q_d = new_q;

        if (A_d.nonZeros() > 0)
        {
            b_d = b_d - A_d * s_vec;
            A_d = A_d * m_vec.asDiagonal();
        }

        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            lb_d(i) = 0.0;
            ub_d(i) = 1.0;
        }
    }

    out.P_rm = (P_d.nonZeros() > 0)
                  ? Eigen::SparseMatrix<double, Eigen::RowMajor>(P_d)
                  : Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n);
    out.A_rm = (A_d.nonZeros() > 0)
                  ? Eigen::SparseMatrix<double, Eigen::RowMajor>(A_d)
                  : Eigen::SparseMatrix<double, Eigen::RowMajor>(static_cast<int>(b_d.size()), n);
    out.q_d = std::move(q_d);
    out.b_d = std::move(b_d);
    out.lb_d = std::move(lb_d);
    out.ub_d = std::move(ub_d);
    out.vtype.assign(n, 'C');
    for (int i = bin_start; i < bin_start + bin_count; ++i)
    {
        out.vtype[i] = 'B';
    }

    return out;
}

} // anonymous namespace

bool gurobi_available()
{
    return GurobiApi::instance().is_available();
}

const std::string& gurobi_unavailable_reason()
{
    return GurobiApi::instance().unavailable_reason();
}

OptSolution solve_qp_gurobi(const Eigen::SparseMatrix<zono_float>& P,
                            const Eigen::Vector<zono_float, -1>& q,
                            zono_float c,
                            const Eigen::SparseMatrix<zono_float>& A,
                            const Eigen::Vector<zono_float, -1>& b,
                            const Eigen::Vector<zono_float, -1>& xi_lb,
                            const Eigen::Vector<zono_float, -1>& xi_ub,
                            const GurobiSettings& settings)
{
    const int n = static_cast<int>(q.size());

    Eigen::SparseMatrix<double, Eigen::RowMajor> P_rm = to_double_rowmajor(P);
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_rm = to_double_rowmajor(A);
    Eigen::VectorXd q_d  = to_double_vec(q);
    Eigen::VectorXd b_d  = to_double_vec(b);
    Eigen::VectorXd lb_d = to_double_vec(xi_lb);
    Eigen::VectorXd ub_d = to_double_vec(xi_ub);

    std::vector<char> vtype(n, 'C');

    PreparedModel pm = prepare_model(P_rm, q_d, A_rm, b_d, lb_d, ub_d, vtype, settings);
    GurobiResult gr = optimize_single(pm);

    Eigen::Vector<zono_float, -1> z = Eigen::Vector<zono_float, -1>::Zero(n);
    if (gr.got_solution)
    {
        z = gr.y.cast<zono_float>();
    }
    return build_opt_solution(gr, z, c);
}

OptSolution solve_miqp_gurobi(const Eigen::SparseMatrix<zono_float>& P,
                              const Eigen::Vector<zono_float, -1>& q,
                              zono_float c,
                              const Eigen::SparseMatrix<zono_float>& A,
                              const Eigen::Vector<zono_float, -1>& b,
                              const Eigen::Vector<zono_float, -1>& xi_lb,
                              const Eigen::Vector<zono_float, -1>& xi_ub,
                              int bin_start, int bin_count,
                              bool zero_one_form,
                              const GurobiSettings& settings)
{
    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, bin_start, bin_count, zero_one_form);

    PreparedModel pm = prepare_model(prep.P_rm, prep.q_d, prep.A_rm, prep.b_d,
                                     prep.lb_d, prep.ub_d, prep.vtype, settings);
    GurobiResult gr = optimize_single(pm);

    const int n = static_cast<int>(q.size());
    Eigen::Vector<zono_float, -1> z = Eigen::Vector<zono_float, -1>::Zero(n);
    if (gr.got_solution)
    {
        Eigen::VectorXd xi = prep.recover_xi(gr.y);
        z = xi.cast<zono_float>();
    }
    return build_opt_solution(gr, z, c + static_cast<zono_float>(prep.const_shift));
}

std::vector<OptSolution> solve_miqp_gurobi_multisol(const Eigen::SparseMatrix<zono_float>& P,
                                                    const Eigen::Vector<zono_float, -1>& q,
                                                    zono_float c,
                                                    const Eigen::SparseMatrix<zono_float>& A,
                                                    const Eigen::Vector<zono_float, -1>& b,
                                                    const Eigen::Vector<zono_float, -1>& xi_lb,
                                                    const Eigen::Vector<zono_float, -1>& xi_ub,
                                                    int bin_start, int bin_count,
                                                    bool zero_one_form,
                                                    int n_sols,
                                                    const GurobiSettings& settings)
{
    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, bin_start, bin_count, zero_one_form);

    PreparedModel pm = prepare_model(prep.P_rm, prep.q_d, prep.A_rm, prep.b_d,
                                     prep.lb_d, prep.ub_d, prep.vtype, settings);
    // If the user explicitly capped PoolSolutions, respect it; otherwise honor n_sols.
    const int capped_n_sols = settings.PoolSolutions
        ? std::min(n_sols, *settings.PoolSolutions)
        : n_sols;
    // Default to PoolSearchMode = 2 ("find n best") unless the user requested otherwise.
    const int pool_search_mode = settings.PoolSearchMode.value_or(2);
    std::vector<GurobiResult> grs = optimize_pool(pm, capped_n_sols, pool_search_mode);

    const int n = static_cast<int>(q.size());
    std::vector<OptSolution> sols;
    sols.reserve(grs.size());
    for (const auto& gr : grs)
    {
        Eigen::Vector<zono_float, -1> z = Eigen::Vector<zono_float, -1>::Zero(n);
        if (gr.got_solution)
        {
            Eigen::VectorXd xi = prep.recover_xi(gr.y);
            z = xi.cast<zono_float>();
        }
        sols.push_back(build_opt_solution(gr, z, c + static_cast<zono_float>(prep.const_shift)));
    }
    return sols;
}

} // namespace detail
} // namespace ZonoOpt
