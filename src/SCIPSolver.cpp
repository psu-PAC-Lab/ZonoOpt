#include "zonoopt/SCIPSolver.hpp"
#include "zonoopt/SCIPApi.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
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

    // Solver-native metadata fetched after each solve for SCIPSolverResults.
    long long node_count = 0;
    double mip_gap = 0.0;
    double dual_bound = 0.0;
    int n_sols_pool = 0;
};

/**
 * RAII wrapper for a fully-built SCIP problem (env + variables + constraints).
 *
 * Member declaration order is intentional: scip_env is declared FIRST so it is
 * destroyed LAST. SCIPreleaseVar / SCIPreleaseCons in the ScopedVar/ScopedCons
 * destructors require the SCIP env to still be valid; with this layout the
 * conss/vars destructors run before scip_env's.
 */
struct ScipProblem
{
    SCIPApi::Scip scip_env;
    std::vector<ScopedVar> vars;       // original n variables, then optionally t_epi
    std::vector<ScopedCons> conss;     // linear equalities, then optionally the quadratic epigraph
    SCIPApi::VarPtr t_var = nullptr;   // auxiliary epigraph variable (non-owning; ScopedVar in vars owns it)
    bool has_quad = false;
    int n = 0;

    SCIPApi::ScipPtr scip() const { return scip_env.get(); }
};

/**
 * Build the SCIP problem (env + vars + linear cons + quadratic epigraph cons if needed)
 * but do not solve. Callers can solve, read solutions, and optionally add more constraints
 * (e.g., no-good cuts) before resolving.
 */
ScipProblem build_scip_problem(const Eigen::SparseMatrix<double>& P_d,
                               const Eigen::VectorXd& q_d,
                               const Eigen::SparseMatrix<double>& A_d,
                               const Eigen::VectorXd& b_d,
                               const Eigen::VectorXd& lb_d,
                               const Eigen::VectorXd& ub_d,
                               const std::vector<char>& vtype,
                               const SCIPSettings& settings)
{
    SCIPApi& api = SCIPApi::instance();
    if (!api.is_available())
        throw std::runtime_error("SCIP API is not available.");

    ScipProblem prob;
    prob.scip_env = api.create_scip();
    prob.n = static_cast<int>(q_d.size());
    prob.has_quad = P_d.nonZeros() > 0;

    SCIPApi::ScipPtr scip = prob.scip();

    // Default to silent unless the user opted in via VerbLevel.
    if (!settings.VerbLevel.has_value() && api.SCIPsetMessagehdlrQuiet)
    {
        api.SCIPsetMessagehdlrQuiet(scip, 1u);
    }

    // Tell SCIP to assume nonlinear constraints are convex. ZonoOpt's QP/MIQP
    // problems always have PSD P (objectives are quadratic distances like ||G xi - x||^2),
    // and the epigraph reformulation we use below produces a convex quadratic constraint.
    // Without this, SCIP's curvature detection on a sparse 35×35 PSD form falls back to
    // generic spatial branching and the solve runs orders of magnitude slower.
    // Applied before apply_scip_settings so a user who knows better can override via the
    // bool_params escape hatch.
    scip_check(api.SCIPsetBoolParam(scip, "constraints/nonlinear/assumeconvex", 1u),
               "set constraints/nonlinear/assumeconvex");

    // For pure (continuous) QPs with a quadratic objective, SCIP's outer-approximation
    // cutting plane separator can fail to close the gap when the optimal objective is
    // near zero: dual stays at 0 (the epigraph variable's lower bound) while primal
    // hovers at some small positive value, leaving the relative gap infinite indefinitely.
    // SCIP's sub-NLP heuristic finds the actual optimum on iteration 1, so terminate
    // after the first feasible solution. Skip this for MIP/MIQP (still need branch-and-bound
    // to enumerate integer assignments) and for pure LP (no quadratic constraint, gap closes
    // immediately at the root, and the first feasible is generally NOT the LP optimum).
    const bool has_binary = std::any_of(vtype.begin(), vtype.end(),
                                        [](char t) { return t == 'B'; });
    if (prob.has_quad && !has_binary)
    {
        scip_check(api.SCIPsetIntParam(scip, "limits/solutions", 1),
                   "set limits/solutions");
    }

    apply_scip_settings(api, scip, settings);

    scip_check(api.SCIPcreateProbBasic(scip, "zonoopt_scip"), "createProbBasic");

    // ---- Variables -------------------------------------------------------
    prob.vars.reserve(prob.n + (prob.has_quad ? 1 : 0));
    for (int i = 0; i < prob.n; ++i)
    {
        const int kind = (vtype[i] == 'B') ? SCIPApi::SCIP_VARTYPE_BINARY
                                           : SCIPApi::SCIP_VARTYPE_CONTINUOUS;
        SCIPApi::VarPtr v = nullptr;
        const std::string name = "x" + std::to_string(i);
        scip_check(api.SCIPcreateVarBasic(scip, &v, name.c_str(),
                                          lb_d(i), ub_d(i), q_d(i), kind),
                   "createVarBasic");
        scip_check(api.SCIPaddVar(scip, v), "addVar");
        prob.vars.emplace_back(scip, v);
    }

    // Epigraph variable t (coefficient 1 in objective) — only if we have a quadratic part.
    // Lower bound 0 is valid because P is PSD (we set assumeconvex above), hence
    // 0.5 xi^T P xi >= 0, hence t = 0.5 xi^T P xi >= 0 at any feasible point. Without
    // this bound the LP relaxation's dual is -inf and SCIP can't close the gap on
    // near-zero QPs via cuts in any reasonable time.
    if (prob.has_quad)
    {
        const double inf = api.SCIPinfinity(scip);
        scip_check(api.SCIPcreateVarBasic(scip, &prob.t_var, "t_epi", 0.0, inf, 1.0,
                                          SCIPApi::SCIP_VARTYPE_CONTINUOUS),
                   "createVarBasic(t_epi)");
        scip_check(api.SCIPaddVar(scip, prob.t_var), "addVar(t_epi)");
        prob.vars.emplace_back(scip, prob.t_var);
    }

    // ---- Linear equality constraints -------------------------------------
    if (A_d.nonZeros() > 0)
    {
        Eigen::SparseMatrix<double, Eigen::RowMajor> A_rm(A_d);
        for (int row = 0; row < A_rm.outerSize(); ++row)
        {
            std::vector<SCIPApi::VarPtr> linvars;
            std::vector<double> lincoefs;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A_rm, row); it; ++it)
            {
                linvars.push_back(prob.vars[it.col()].var);
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
            prob.conss.emplace_back(scip, c);
        }
    }

    // ---- Quadratic epigraph constraint: 0.5 * xi^T P xi - t <= 0 ---------
    if (prob.has_quad)
    {
        if (!api.SCIPcreateConsBasicQuadraticNonlinear)
        {
            throw std::runtime_error(
                "SCIPcreateConsBasicQuadraticNonlinear is unavailable in this SCIP build; "
                "SCIP 8 or newer is required for QP/MIQP support via ZonoOpt.");
        }

        SCIPApi::VarPtr lvars[1] = { prob.t_var };
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
                if (it.col() > it.row()) continue;
                const double v = (it.row() == it.col()) ? 0.5 * it.value() : it.value();
                qv1.push_back(prob.vars[it.row()].var);
                qv2.push_back(prob.vars[it.col()].var);
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
        prob.conss.emplace_back(scip, qc_cons);
    }

    return prob;
}

/**
 * Solve `prob` once and extract the best solution into a SolveResult.
 * On infeasibility the SolveResult has infeasible=true and got_solution=false.
 */
SolveResult solve_once(ScipProblem& prob)
{
    SCIPApi& api = SCIPApi::instance();
    SCIPApi::ScipPtr scip = prob.scip();

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
        res.y = Eigen::VectorXd::Zero(prob.n);
        for (int i = 0; i < prob.n; ++i)
            res.y(i) = api.SCIPgetSolVal(scip, best, prob.vars[i].var);
        res.objective = api.SCIPgetSolOrigObj ? api.SCIPgetSolOrigObj(scip, best) : 0.0;
        res.got_solution = true;
    }

    // Solver-native metadata for SCIPSolverResults (best-effort; symbols may be absent on older builds).
    if (api.SCIPgetNTotalNodes) res.node_count  = api.SCIPgetNTotalNodes(scip);
    if (api.SCIPgetGap)         res.mip_gap     = api.SCIPgetGap(scip);
    if (api.SCIPgetDualbound)   res.dual_bound  = api.SCIPgetDualbound(scip);
    if (api.SCIPgetNSols)       res.n_sols_pool = api.SCIPgetNSols(scip);

    return res;
}

/**
 * Backward-compatible single-shot solve used by solve_qp_scip / solve_miqp_scip.
 * (Kept as a thin wrapper for clarity at the call sites.)
 */
SolveResult run_scip(const Eigen::SparseMatrix<double>& P_d,
                     const Eigen::VectorXd& q_d,
                     const Eigen::SparseMatrix<double>& A_d,
                     const Eigen::VectorXd& b_d,
                     const Eigen::VectorXd& lb_d,
                     const Eigen::VectorXd& ub_d,
                     const std::vector<char>& vtype,
                     int /*n_sols, unused*/,
                     const SCIPSettings& settings)
{
    ScipProblem prob = build_scip_problem(P_d, q_d, A_d, b_d, lb_d, ub_d, vtype, settings);
    return solve_once(prob);
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

    // Attach SCIP-native solution metadata so callers can read raw status, node count, etc.
    auto sres = std::make_shared<SCIPSolverResults>();
    sres->status      = sr.status;
    sres->node_count  = sr.node_count;
    sres->mip_gap     = sr.mip_gap;
    sres->dual_bound  = sr.dual_bound;
    sres->n_sols_pool = sr.n_sols_pool;
    sol.external_results = sres;

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

    SCIPApi& api = SCIPApi::instance();
    if (!api.is_available())
        throw std::runtime_error("SCIP API is not available.");

    MiqpPrep prep = prep_miqp(P, q, A, b, xi_lb, xi_ub, bin_start, bin_count, zero_one_form);
    const int n = static_cast<int>(q.size());
    const zono_float offset = c + static_cast<zono_float>(prep.const_shift);

    // With no binary variables there's nothing to enumerate; return the QP optimum.
    if (bin_count == 0)
    {
        return { solve_qp_scip(P, q, c, A, b, xi_lb, xi_ub, settings) };
    }

    if (!api.SCIPfreeTransform)
        throw std::runtime_error("SCIPfreeTransform is unavailable; cannot enumerate multiple feasible binary assignments.");

    // Cap requested count at 2^bin_count (the maximum number of distinct binary assignments).
    {
        const int max_combinations =
            (bin_count < 31) ? (1 << bin_count) : (std::numeric_limits<int>::max)();
        n_sols = (std::min)(n_sols, max_combinations);
    }

    // Enumerate feasible binary assignments via repeated solve + no-good cuts.
    //
    // The previous implementation used SCIP's count constraint handler, but that handler
    // populates its variable list via SCIPgetVarsData during solving, which returns
    // *transformed* variable pointers. Our prob.vars array holds *original* variable
    // pointers (from SCIPcreateVarBasic + SCIPaddVar). The pointer-based lookup always
    // failed, leaving all binary values at zero — infeasible for union sets whose
    // indicator variables must sum to one.
    //
    // The no-good cut approach avoids the issue entirely: SCIPgetSolVal handles
    // original/transformed translation internally, so sr.y always carries the correct
    // binary values read from the current solution.

    // Feasibility model: zero objective (any feasible binary assignment is acceptable).
    Eigen::SparseMatrix<double> P_empty(prep.P_d.rows(), prep.P_d.cols());
    Eigen::VectorXd q_zero = Eigen::VectorXd::Zero(prep.q_d.size());

    ScipProblem prob = build_scip_problem(P_empty, q_zero, prep.A_d, prep.b_d,
                                          prep.lb_d, prep.ub_d, prep.vtype, settings);
    SCIPApi::ScipPtr scip = prob.scip();
    const double inf = api.SCIPinfinity(scip);

    std::vector<OptSolution> sols;
    // Do not reserve based on n_sols — it may be INT_MAX for "find all leaves".

    const auto t_global_start = std::chrono::steady_clock::now();
    bool converged = false;

    for (int k = 0; k < n_sols; ++k)
    {
        SolveResult sr = solve_once(prob);

        if (sr.infeasible)
        {
            // All feasible binary assignments have been enumerated.
            converged = true;
            break;
        }
        if (!sr.got_solution)
            break;  // time limit or other early termination without a solution

        // Read binary values from sr.y (populated by SCIPgetSolVal in solve_once, which
        // handles original↔transformed variable translation internally). Round to nearest
        // integer to guard against small floating-point errors near 0 or 1.
        Eigen::VectorXd y = Eigen::VectorXd::Zero(prep.q_d.size());
        for (int i = bin_start; i < bin_start + bin_count; ++i)
            y(i) = std::round(sr.y(i));

        OptSolution sol;
        sol.z = prep.recover_xi(y).cast<zono_float>();
        sol.J = offset;
        sol.converged = true;  // overwritten on the last entry after the loop
        sol.infeasible = false;
        sol.iter = 0;
        sol.run_time = sr.runtime;
        {
            auto sres = std::make_shared<SCIPSolverResults>();
            sres->status      = sr.status;
            sres->node_count  = sr.node_count;
            sres->mip_gap     = 0.0;
            sres->dual_bound  = 0.0;
            sres->n_sols_pool = static_cast<int>(sols.size()) + 1;
            sol.external_results = sres;
        }
        sols.push_back(std::move(sol));

        if (static_cast<int>(sols.size()) == n_sols)
            break;

        // Free the transformed problem so a new constraint can be added to the original.
        scip_check(api.SCIPfreeTransform(scip), "SCIPfreeTransform");

        // No-good cut: at least one binary variable must differ from the current assignment.
        // In y-space (all binaries in {0,1}):
        //   sum_{i: y_i=1} (1-x_i) + sum_{i: y_i=0} x_i >= 1
        //   => coef_i = -1 if y_i=1 else +1,  lhs = 1 - count(y_i=1)
        std::vector<SCIPApi::VarPtr> ng_vars;
        std::vector<double> ng_coefs;
        ng_vars.reserve(static_cast<size_t>(bin_count));
        ng_coefs.reserve(static_cast<size_t>(bin_count));
        double ng_lhs = 1.0;
        for (int i = bin_start; i < bin_start + bin_count; ++i)
        {
            ng_vars.push_back(prob.vars[i].var);
            if (y(i) > 0.5)
            {
                ng_coefs.push_back(-1.0);
                ng_lhs -= 1.0;
            }
            else
            {
                ng_coefs.push_back(1.0);
            }
        }
        SCIPApi::ConsPtr ng_cons = nullptr;
        const std::string ng_name = "nogood_" + std::to_string(k);
        scip_check(api.SCIPcreateConsBasicLinear(scip, &ng_cons, ng_name.c_str(),
                                                  static_cast<int>(ng_vars.size()),
                                                  ng_vars.data(), ng_coefs.data(),
                                                  ng_lhs, inf),
                   "createConsBasicLinear(nogood)");
        scip_check(api.SCIPaddCons(scip, ng_cons), "addCons(nogood)");
        scip_check(api.SCIPreleaseCons(scip, &ng_cons), "releaseCons(nogood)");

        // Update the per-solve time limit so total enumeration respects the budget.
        if (settings.TimeLimit)
        {
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_global_start).count();
            const double remaining = *settings.TimeLimit - elapsed;
            if (remaining <= 0.0)
                break;
            scip_check(api.SCIPsetRealParam(scip, "limits/time", remaining),
                       "set limits/time");
        }
    }

    if (sols.empty())
    {
        // Infeasible on the very first solve — the set has no feasible binary assignment.
        OptSolution s;
        s.z = Eigen::Vector<zono_float, -1>::Zero(n);
        s.infeasible = true;
        s.converged = false;
        s.J = std::numeric_limits<zono_float>::infinity();
        s.run_time = 0.0;
        s.iter = 0;
        {
            auto sres = std::make_shared<SCIPSolverResults>();
            sres->status      = SCIPApi::SCIP_STATUS_INFEASIBLE;
            sres->node_count  = 0;
            sres->mip_gap     = 0.0;
            sres->dual_bound  = 0.0;
            sres->n_sols_pool = 0;
            s.external_results = sres;
        }
        return { std::move(s) };
    }

    // Mark the last solution with the overall convergence status.
    sols.back().converged = converged;
    return sols;
}

} // namespace detail
} // namespace ZonoOpt
