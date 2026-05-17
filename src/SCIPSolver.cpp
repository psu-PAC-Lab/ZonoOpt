#include "zonoopt/SCIPSolver.hpp"
#include "zonoopt/SCIPApi.hpp"

#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ZonoOpt
{
namespace detail
{

namespace
{

// ---- Small helpers ---------------------------------------------------------

inline void scip_check(int rc, const char* what)
{
    if (rc != SCIPApi::SCIP_OKAY)
    {
        throw std::runtime_error(std::string("SCIP ") + what + " failed (rc=" + std::to_string(rc) + ").");
    }
}

template <typename T>
Eigen::VectorXd to_double_vec(const Eigen::Vector<T, -1>& v) { return v.template cast<double>(); }

/// RAII wrapper for SCIP_VAR* — calls SCIPreleaseVar on destruction.
struct ScopedVar
{
    SCIPApi::ScipPtr scip = nullptr;
    SCIPApi::VarPtr  var  = nullptr;
    ScopedVar() = default;
    ScopedVar(SCIPApi::ScipPtr s, SCIPApi::VarPtr v) : scip(s), var(v) {}
    ScopedVar(const ScopedVar&) = delete;
    ScopedVar& operator=(const ScopedVar&) = delete;
    ScopedVar(ScopedVar&& other) noexcept : scip(other.scip), var(other.var) { other.var = nullptr; }
    ~ScopedVar()
    {
        if (var && scip)
        {
            auto& api = SCIPApi::instance();
            if (api.SCIPreleaseVar) api.SCIPreleaseVar(scip, &var);
        }
    }
};

/// RAII wrapper for SCIP_CONS* — calls SCIPreleaseCons on destruction.
struct ScopedCons
{
    SCIPApi::ScipPtr scip = nullptr;
    SCIPApi::ConsPtr cons = nullptr;
    ScopedCons() = default;
    ScopedCons(SCIPApi::ScipPtr s, SCIPApi::ConsPtr c) : scip(s), cons(c) {}
    ScopedCons(const ScopedCons&) = delete;
    ScopedCons& operator=(const ScopedCons&) = delete;
    ScopedCons(ScopedCons&& other) noexcept : scip(other.scip), cons(other.cons) { other.cons = nullptr; }
    ~ScopedCons()
    {
        if (cons && scip)
        {
            auto& api = SCIPApi::instance();
            if (api.SCIPreleaseCons) api.SCIPreleaseCons(scip, &cons);
        }
    }
};

/// Apply all parameters set by the user. Maps onto SCIPset{Bool,Int,Longint,Real,Char,String}Param.
void apply_scip_settings(SCIPApi& api, SCIPApi::ScipPtr scip, const SCIPSettings& s)
{
    auto apply_int = [&](const char* name, const std::optional<int>& v) {
        if (v) scip_check(api.SCIPsetIntParam(scip, name, *v), "set int param");
    };
    auto apply_dbl = [&](const char* name, const std::optional<double>& v) {
        if (v) scip_check(api.SCIPsetRealParam(scip, name, *v), "set real param");
    };

    // Typed fields → SCIP parameter names.
    apply_dbl("limits/time",                       s.TimeLimit);
    apply_dbl("limits/memory",                     s.MemLimit);
    apply_int("limits/solutions",                  s.SolutionLimit);
    if (s.NodeLimit)
    {
        // NodeLimit is a long long in SCIP; expose as int here for ergonomics.
        scip_check(api.SCIPsetLongintParam(scip, "limits/nodes",
                                           static_cast<long long>(*s.NodeLimit)),
                   "set limits/nodes");
    }
    apply_dbl("limits/gap",                        s.MIPGap);
    apply_dbl("limits/absgap",                     s.MIPGapAbs);
    apply_dbl("numerics/feastol",                  s.FeasibilityTol);
    apply_int("parallel/maxnthreads",              s.Threads);
    apply_int("randomization/randomseedshift",     s.Seed);
    apply_int("display/verblevel",                 s.VerbLevel);

    // Escape-hatch maps.
    for (const auto& kv : s.bool_params)
        scip_check(api.SCIPsetBoolParam(scip, kv.first.c_str(), kv.second ? 1u : 0u), "set bool param");
    for (const auto& kv : s.int_params)
        scip_check(api.SCIPsetIntParam(scip, kv.first.c_str(), kv.second), "set int param");
    for (const auto& kv : s.longint_params)
        scip_check(api.SCIPsetLongintParam(scip, kv.first.c_str(), kv.second), "set longint param");
    for (const auto& kv : s.real_params)
        scip_check(api.SCIPsetRealParam(scip, kv.first.c_str(), kv.second), "set real param");
    for (const auto& kv : s.char_params)
        scip_check(api.SCIPsetCharParam(scip, kv.first.c_str(), kv.second), "set char param");
    for (const auto& kv : s.str_params)
        scip_check(api.SCIPsetStringParam(scip, kv.first.c_str(), kv.second.c_str()), "set string param");
}

// ---- {-1,1} → {0,1} substitution for MIQP binaries (mirrors GurobiSolver) -----
//
// When binary variables xi_i are in {-1, 1}, SCIP (like Gurobi) only supports
// {0, 1} binaries directly. We substitute xi_i = 2*y_i - 1 (only on binary indices)
// before building the SCIP model, then undo the substitution on the recovered
// solution.

struct MiqpPrep
{
    // After substitution (in y-space when needs_transform; otherwise unchanged).
    Eigen::SparseMatrix<double> P_d;
    Eigen::SparseMatrix<double> A_d;
    Eigen::VectorXd q_d;
    Eigen::VectorXd b_d;
    Eigen::VectorXd lb_d;
    Eigen::VectorXd ub_d;
    // Variable kind: 'B' for binary in {0,1}, 'C' continuous, with the original
    // continuous lb/ub for non-binary indices.
    std::vector<char> vtype;

    double const_shift = 0.0;
    bool needs_transform = false;
    int bin_start = 0;
    int bin_count = 0;

    // Convert a solution from y-space back to xi-space.
    Eigen::VectorXd recover_xi(const Eigen::VectorXd& y) const
    {
        if (!needs_transform) return y;
        Eigen::VectorXd xi = y;
        for (int i = bin_start; i < bin_start + bin_count; ++i)
            xi(i) = 2.0 * y(i) - 1.0;
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

    MiqpPrep out;
    out.bin_start = bin_start;
    out.bin_count = bin_count;
    out.needs_transform = (bin_count > 0) && !zero_one_form;

    out.P_d = (P.nonZeros() > 0)
                ? Eigen::SparseMatrix<double>(P.template cast<double>())
                : Eigen::SparseMatrix<double>(n, n);
    out.A_d = (A.nonZeros() > 0)
                ? Eigen::SparseMatrix<double>(A.template cast<double>())
                : Eigen::SparseMatrix<double>(static_cast<int>(b.size()), n);
    out.q_d  = to_double_vec(q);
    out.b_d  = to_double_vec(b);
    out.lb_d = to_double_vec(xi_lb);
    out.ub_d = to_double_vec(xi_ub);

    if (out.needs_transform)
    {
        Eigen::VectorXd m_vec = Eigen::VectorXd::Ones(n);
        Eigen::VectorXd s_vec = Eigen::VectorXd::Zero(n);
        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            m_vec(i) = 2.0;
            s_vec(i) = -1.0;
        }

        Eigen::VectorXd P_s = (out.P_d.nonZeros() > 0)
                                ? Eigen::VectorXd(out.P_d * s_vec)
                                : Eigen::VectorXd::Zero(n);
        Eigen::VectorXd new_q = m_vec.asDiagonal() * (P_s + out.q_d);
        out.const_shift = 0.5 * s_vec.dot(P_s) + out.q_d.dot(s_vec);

        if (out.P_d.nonZeros() > 0)
            out.P_d = m_vec.asDiagonal() * out.P_d * m_vec.asDiagonal();
        out.q_d = new_q;

        if (out.A_d.nonZeros() > 0)
        {
            out.b_d = out.b_d - out.A_d * s_vec;
            out.A_d = out.A_d * m_vec.asDiagonal();
        }

        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            out.lb_d(i) = 0.0;
            out.ub_d(i) = 1.0;
        }
    }

    out.vtype.assign(n, 'C');
    for (int i = bin_start; i < bin_start + bin_count; ++i)
        out.vtype[i] = 'B';

    return out;
}

// ---- Model build / solve / extract -----------------------------------------

struct SolveResult
{
    bool got_solution = false;
    bool infeasible = false;
    bool optimal = false;
    int status = 0;
    double objective = 0.0;
    double runtime = 0.0;
    Eigen::VectorXd y;            // primary solution (best)
    std::vector<Eigen::VectorXd> extra_solutions;  // additional solutions (pool); each is in y-space
    std::vector<double>          extra_objectives;
};

/**
 * Build a SCIP model from the given (already y-space, if applicable) problem data,
 * solve it, and extract up to `n_sols` feasible solutions.
 *
 * On `n_sols > 1`, "limits/maxsol" is set so SCIP retains that many solutions during
 * search, and after solve we copy them out via SCIPgetSols / SCIPgetNSols.
 */
SolveResult run_scip(const Eigen::SparseMatrix<double>& P_d,
                     const Eigen::VectorXd& q_d,
                     const Eigen::SparseMatrix<double>& A_d,
                     const Eigen::VectorXd& b_d,
                     const Eigen::VectorXd& lb_d,
                     const Eigen::VectorXd& ub_d,
                     const std::vector<char>& vtype,
                     int n_sols,
                     const SCIPSettings& settings)
{
    SCIPApi& api = SCIPApi::instance();
    if (!api.is_available())
        throw std::runtime_error("SCIP API is not available.");

    SCIPApi::Scip scip_env = api.create_scip();
    SCIPApi::ScipPtr scip = scip_env.get();

    // Default to silent unless the user opted in via VerbLevel.
    if (!settings.VerbLevel.has_value() && api.SCIPsetMessagehdlrQuiet)
    {
        api.SCIPsetMessagehdlrQuiet(scip, 1u);
    }

    // Multisol storage capacity. SCIP's solution pool holds up to "limits/maxsol"
    // feasible solutions found during search.
    if (n_sols > 1)
    {
        scip_check(api.SCIPsetIntParam(scip, "limits/maxsol", std::max(n_sols, 1)),
                   "set limits/maxsol");
    }

    apply_scip_settings(api, scip, settings);

    scip_check(api.SCIPcreateProbBasic(scip, "zonoopt_scip"), "createProbBasic");

    const int n = static_cast<int>(q_d.size());
    const bool has_quad = P_d.nonZeros() > 0;

    // ---- Variables -------------------------------------------------------
    std::vector<ScopedVar> vars;
    vars.reserve(n + (has_quad ? 1 : 0));
    for (int i = 0; i < n; ++i)
    {
        const int kind = (vtype[i] == 'B') ? SCIPApi::SCIP_VARTYPE_BINARY
                                           : SCIPApi::SCIP_VARTYPE_CONTINUOUS;
        SCIPApi::VarPtr v = nullptr;
        const std::string name = "x" + std::to_string(i);
        scip_check(api.SCIPcreateVarBasic(scip, &v, name.c_str(),
                                          lb_d(i), ub_d(i), q_d(i), kind),
                   "createVarBasic");
        scip_check(api.SCIPaddVar(scip, v), "addVar");
        vars.emplace_back(scip, v);
    }

    // Epigraph variable t (coefficient 1 in objective) — only if we have a quadratic part.
    SCIPApi::VarPtr t_var = nullptr;
    if (has_quad)
    {
        const double inf = api.SCIPinfinity(scip);
        scip_check(api.SCIPcreateVarBasic(scip, &t_var, "t_epi", -inf, inf, 1.0,
                                          SCIPApi::SCIP_VARTYPE_CONTINUOUS),
                   "createVarBasic(t_epi)");
        scip_check(api.SCIPaddVar(scip, t_var), "addVar(t_epi)");
        vars.emplace_back(scip, t_var);
    }

    // ---- Linear equality constraints -------------------------------------
    std::vector<ScopedCons> conss;
    if (A_d.nonZeros() > 0)
    {
        // A is column-major; iterate rows by walking the row-major copy.
        Eigen::SparseMatrix<double, Eigen::RowMajor> A_rm(A_d);
        for (int row = 0; row < A_rm.outerSize(); ++row)
        {
            std::vector<SCIPApi::VarPtr> linvars;
            std::vector<double> lincoefs;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_rm, row); it; ++it)
            {
                linvars.push_back(vars[it.col()].var);
                lincoefs.push_back(it.value());
            }
            if (linvars.empty()) continue;
            SCIPApi::ConsPtr c = nullptr;
            const std::string name = "eq" + std::to_string(row);
            const double rhs = b_d(row);
            scip_check(api.SCIPcreateConsBasicLinear(scip, &c, name.c_str(),
                                                     static_cast<int>(linvars.size()),
                                                     linvars.data(), lincoefs.data(),
                                                     rhs, rhs),
                       "createConsBasicLinear");
            scip_check(api.SCIPaddCons(scip, c), "addCons");
            conss.emplace_back(scip, c);
        }
    }

    // ---- Quadratic constraint for the objective epigraph -----------------
    //
    // 0.5 * xi^T P xi - t <= 0
    //
    // Build the term list from P (symmetric, lower-triangular contribution scaled
    // by 0.5 on the diagonal, full off-diagonal counted once).
    if (has_quad)
    {
        if (!api.SCIPcreateConsBasicQuadraticNonlinear)
        {
            throw std::runtime_error(
                "SCIPcreateConsBasicQuadraticNonlinear is unavailable in this SCIP build; "
                "SCIP 8 or newer is required for QP/MIQP support via ZonoOpt.");
        }

        // Linear part: -1 * t
        SCIPApi::VarPtr lvars[1] = { t_var };
        double          lcoefs[1] = { -1.0 };

        std::vector<SCIPApi::VarPtr> qv1, qv2;
        std::vector<double> qc;
        qv1.reserve(P_d.nonZeros());
        qv2.reserve(P_d.nonZeros());
        qc.reserve(P_d.nonZeros());

        for (int k = 0; k < P_d.outerSize(); ++k)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(P_d, k); it; ++it)
            {
                if (it.col() > it.row()) continue; // lower-triangular only
                const double v = (it.row() == it.col()) ? 0.5 * it.value() : it.value();
                qv1.push_back(vars[it.row()].var);
                qv2.push_back(vars[it.col()].var);
                qc.push_back(v);
            }
        }

        const double inf = api.SCIPinfinity(scip);
        SCIPApi::ConsPtr qc_cons = nullptr;
        scip_check(api.SCIPcreateConsBasicQuadraticNonlinear(
                       scip, &qc_cons, "obj_epigraph",
                       1, lvars, lcoefs,
                       static_cast<int>(qv1.size()), qv1.data(), qv2.data(), qc.data(),
                       -inf, 0.0),
                   "createConsBasicQuadraticNonlinear");
        scip_check(api.SCIPaddCons(scip, qc_cons), "addCons(quadratic)");
        conss.emplace_back(scip, qc_cons);
    }

    // ---- Solve -----------------------------------------------------------
    auto t_start = std::chrono::steady_clock::now();
    scip_check(api.SCIPsolve(scip), "SCIPsolve");
    auto t_end = std::chrono::steady_clock::now();

    SolveResult res;
    res.runtime = std::chrono::duration<double>(t_end - t_start).count();
    res.status  = api.SCIPgetStatus(scip);
    res.optimal = (res.status == api.SCIP_STATUS_OPTIMAL);
    res.infeasible = (res.status == api.SCIP_STATUS_INFEASIBLE
                   || res.status == api.SCIP_STATUS_INFORUNBD);

    SCIPApi::SolPtr best = api.SCIPgetBestSol(scip);
    if (best != nullptr && !res.infeasible)
    {
        res.y = Eigen::VectorXd::Zero(n);
        for (int i = 0; i < n; ++i)
            res.y(i) = api.SCIPgetSolVal(scip, best, vars[i].var);
        res.objective = api.SCIPgetSolOrigObj ? api.SCIPgetSolOrigObj(scip, best) : 0.0;
        res.got_solution = true;
    }

    if (n_sols > 1 && api.SCIPgetSols && api.SCIPgetNSols)
    {
        const int n_in_pool = api.SCIPgetNSols(scip);
        SCIPApi::SolPtr* pool = api.SCIPgetSols(scip);
        const int take = std::min(n_sols, n_in_pool);
        // The "best" pool entry is index 0; we already captured it above. Include
        // the remainder as extras.
        for (int k = 1; k < take; ++k)
        {
            Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
            for (int i = 0; i < n; ++i)
                y(i) = api.SCIPgetSolVal(scip, pool[k], vars[i].var);
            res.extra_solutions.push_back(std::move(y));
            res.extra_objectives.push_back(api.SCIPgetSolOrigObj ? api.SCIPgetSolOrigObj(scip, pool[k]) : 0.0);
        }
    }

    return res;
}

OptSolution build_opt_solution(const SolveResult& sr, const Eigen::Vector<zono_float, -1>& z,
                               zono_float constant_offset)
{
    OptSolution sol;
    sol.z = z;
    sol.run_time = sr.runtime;
    sol.iter = 0;

    if (sr.infeasible)
    {
        sol.infeasible = true;
        sol.converged = false;
        sol.J = std::numeric_limits<zono_float>::infinity();
    }
    else if (sr.got_solution)
    {
        sol.infeasible = false;
        sol.converged = sr.optimal;
        sol.J = static_cast<zono_float>(sr.objective) + constant_offset;
    }
    else
    {
        sol.infeasible = false;
        sol.converged = false;
        sol.J = -std::numeric_limits<zono_float>::infinity();
    }
    return sol;
}

} // anonymous namespace

bool scip_available()
{
    return SCIPApi::instance().is_available();
}

const std::string& scip_unavailable_reason()
{
    return SCIPApi::instance().unavailable_reason();
}

// ---------- Public API -------------------------------------------------------

OptSolution solve_qp_scip(const Eigen::SparseMatrix<zono_float>& P,
                          const Eigen::Vector<zono_float, -1>& q,
                          zono_float c,
                          const Eigen::SparseMatrix<zono_float>& A,
                          const Eigen::Vector<zono_float, -1>& b,
                          const Eigen::Vector<zono_float, -1>& xi_lb,
                          const Eigen::Vector<zono_float, -1>& xi_ub,
                          const SCIPSettings& settings)
{
    // Reuse the MIQP prep with 0 binary variables — it just casts to double.
    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, /*bin_start=*/0, /*bin_count=*/0,
                              /*zero_one_form=*/true);

    SolveResult sr = run_scip(prep.P_d, prep.q_d, prep.A_d, prep.b_d,
                              prep.lb_d, prep.ub_d, prep.vtype, /*n_sols=*/1, settings);

    const int n = static_cast<int>(q.size());
    Eigen::Vector<zono_float, -1> z = Eigen::Vector<zono_float, -1>::Zero(n);
    if (sr.got_solution) z = sr.y.cast<zono_float>();
    return build_opt_solution(sr, z, c + static_cast<zono_float>(prep.const_shift));
}

OptSolution solve_miqp_scip(const Eigen::SparseMatrix<zono_float>& P,
                            const Eigen::Vector<zono_float, -1>& q,
                            zono_float c,
                            const Eigen::SparseMatrix<zono_float>& A,
                            const Eigen::Vector<zono_float, -1>& b,
                            const Eigen::Vector<zono_float, -1>& xi_lb,
                            const Eigen::Vector<zono_float, -1>& xi_ub,
                            int bin_start, int bin_count,
                            bool zero_one_form,
                            const SCIPSettings& settings)
{
    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, bin_start, bin_count, zero_one_form);
    SolveResult sr = run_scip(prep.P_d, prep.q_d, prep.A_d, prep.b_d,
                              prep.lb_d, prep.ub_d, prep.vtype, /*n_sols=*/1, settings);

    const int n = static_cast<int>(q.size());
    Eigen::Vector<zono_float, -1> z = Eigen::Vector<zono_float, -1>::Zero(n);
    if (sr.got_solution) z = prep.recover_xi(sr.y).cast<zono_float>();
    return build_opt_solution(sr, z, c + static_cast<zono_float>(prep.const_shift));
}

std::vector<OptSolution> solve_miqp_scip_multisol(const Eigen::SparseMatrix<zono_float>& P,
                                                  const Eigen::Vector<zono_float, -1>& q,
                                                  zono_float c,
                                                  const Eigen::SparseMatrix<zono_float>& A,
                                                  const Eigen::Vector<zono_float, -1>& b,
                                                  const Eigen::Vector<zono_float, -1>& xi_lb,
                                                  const Eigen::Vector<zono_float, -1>& xi_ub,
                                                  int bin_start, int bin_count,
                                                  bool zero_one_form,
                                                  int n_sols,
                                                  const SCIPSettings& settings)
{
    if (n_sols < 1) n_sols = 1;

    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, bin_start, bin_count, zero_one_form);
    SolveResult sr = run_scip(prep.P_d, prep.q_d, prep.A_d, prep.b_d,
                              prep.lb_d, prep.ub_d, prep.vtype, n_sols, settings);

    const int n = static_cast<int>(q.size());
    const zono_float offset = c + static_cast<zono_float>(prep.const_shift);

    std::vector<OptSolution> sols;

    if (sr.infeasible)
    {
        OptSolution s;
        s.z = Eigen::Vector<zono_float, -1>::Zero(n);
        s.infeasible = true;
        s.converged = false;
        s.J = std::numeric_limits<zono_float>::infinity();
        s.run_time = sr.runtime;
        sols.push_back(std::move(s));
        return sols;
    }

    // Best solution first.
    if (sr.got_solution)
    {
        Eigen::Vector<zono_float, -1> z = prep.recover_xi(sr.y).cast<zono_float>();
        sols.push_back(build_opt_solution(sr, z, offset));
    }
    // Pool extras.
    for (size_t k = 0; k < sr.extra_solutions.size(); ++k)
    {
        SolveResult sr_k = sr;            // copy meta (runtime, status, etc.)
        sr_k.y = sr.extra_solutions[k];
        sr_k.objective = sr.extra_objectives[k];
        sr_k.got_solution = true;
        sr_k.optimal = false;             // these are feasible-but-not-proved-optimal
        Eigen::Vector<zono_float, -1> z = prep.recover_xi(sr_k.y).cast<zono_float>();
        sols.push_back(build_opt_solution(sr_k, z, offset));
    }
    return sols;
}

} // namespace detail
} // namespace ZonoOpt
