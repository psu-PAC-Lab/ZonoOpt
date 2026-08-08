import numpy as np
from scipy import sparse
import zonoopt as zono
try:
    import objective_feasibility_pump as ofp
except ImportError as e:
    print("Need to pip install OFP from https://github.com/jrobbins11/objective-feasibility-pump")
    raise e
import json

def random_sparse_matrix(m, n, density, val_range):
    assert(density > 0 and density <= 1)
    assert(m > 0 and n > 0)
    
    max_nnz = m*n
    nnz = int(np.round(max_nnz * density))

    # triplets
    rows = []
    cols = []
    vals = []
    for _ in range(nnz):
        rows.append(np.round(np.random.random()*(m-1)))
        cols.append(np.round(np.random.random()*(n-1)))
        vals.append(np.random.random()*(2*val_range) - val_range)

    # matrix
    return sparse.csc_matrix((vals, (rows, cols)), shape=(m,n))

def random_vector(m, val_range):
    return np.random.rand(m)*val_range - val_range

def random_hz(n, nGc, nGb, nC, density):
    Gc = random_sparse_matrix(n, nGc, density, 1.0)
    Gb = random_sparse_matrix(n, nGb, density, 1.0)
    c = random_vector(n, 1.0)
    Ac = random_sparse_matrix(nC, nGc, density, 1.0)
    Ab = random_sparse_matrix(nC, nGb, density, 1.0)
    b = random_vector(nC, 1.0)

    Z = zono.HybZono(Gc, Gb, c, Ac, Ab, b)
    if not Z.is_0_1_form():
        Z.convert_form()
    Z.remove_redundancy()
    return Z

def admm_fp_solve(Z, P, q, c, settings, seed):
    admm_fp_settings = settings.copy()
    admm_fp_settings.use_interval_contractor = False
    admm_fp_settings.single_threaded_admm_fp = True
    admm_fp_settings.polish = False
    admm_fp_settings.k_max_admm_fp_ph2 = int(1e7) # massive number, just let max time be limiting factor
    admm_fp_settings.enable_restart_admm_fp = True
    admm_fp_settings.enable_perturb_admm_fp = True
    admm_fp_settings.enable_rng_seed = True
    admm_fp_settings.rng_seed = seed

    sol = zono.OptSolution()
    x = Z.optimize_over(P, q, c=c, settings=admm_fp_settings, solution=sol)

    return sol.converged, sol.run_time, sol.z, sol.iter

def admm_solve(Z, P, q, c, settings, seed):
    admm_settings = settings.copy()
    admm_settings.use_interval_contractor = False
    admm_settings.single_threaded_admm_fp = True
    admm_settings.polish = False
    admm_settings.k_max_admm_fp_ph1 = int(1e7) # massive number, just let max time be limiting factor
    admm_settings.enable_restart_admm_fp = False
    admm_settings.enable_perturb_admm_fp = False
    admm_settings.enable_rng_seed = True
    admm_settings.rng_seed = seed

    sol = zono.OptSolution()
    x = Z.optimize_over(P, q, c=c, settings=admm_settings, solution=sol)

    return sol.converged, sol.run_time, sol.z, sol.iter

def ofp_solve(Z, P, q, c, settings, seed):
    ofp_settings = ofp.OFP_Settings()
    ofp_settings.t_max = settings.t_max
    ofp_settings.max_iter = 1000000 # very large number
    ofp_settings.verbose = settings.verbose
    ofp_settings.verbosity_interval = 1
    ofp_settings.max_restarts = 100000 # very large number
    ofp_settings.rng_seed = seed
    ofp_settings.inf_norm_conv = True
    ofp_settings.tol = settings.eps_prim

    q_tilde = Z.get_G().transpose().dot(P*Z.get_c() + q)
    bins = [i for i in range(Z.get_nG() - Z.get_nGb(), Z.get_nG())]

    # solve
    if not Z.is_0_1_form():
        raise ValueError("OFP solver requires 0-1 form hybrid zonotope.")
    solver = ofp.OFP_Solver(q_tilde, Z.get_A(), Z.get_b(), Z.get_b(), np.zeros(Z.get_nG()), np.ones(Z.get_nG()), bins, settings=ofp_settings, b=c)
    success = solver.solve()
    info = solver.get_info()

    return success, info.runtime, solver.get_solution(), info.iter

def gurobi_solve(Z, P, q, c, settings, seed, find_global_opt=False, t_max_global_opt=10.0):
    # solve using ZonoOpt's built-in Gurobi solver support
    gurobi_settings = zono.GurobiSettings()
    if find_global_opt:
        gurobi_settings.MIPGap = 0.0 # must find global optimum
        gurobi_settings.TimeLimit = t_max_global_opt
    else:
        gurobi_settings.MIPGap = 1.0 # any feasible solution is acceptable
        gurobi_settings.TimeLimit = settings.t_max
    gurobi_settings.OutputFlag = 1 if settings.verbose else 0
    gurobi_settings.MIPGapAbs = 0.0
    gurobi_settings.Threads = 1
    gurobi_settings.FeasibilityTol = settings.eps_prim

    sol = zono.OptSolution()
    Z.optimize_over(P, q, c=c, settings=gurobi_settings, solution=sol)

    return sol.converged, sol.run_time, sol.z, np.nan

def scip_solve(Z, P, q, c, settings, seed, find_global_opt=False, t_max_global_opt=10.0):
    # solve using ZonoOpt's built-in SCIP solver support
    scip_settings = zono.SCIPSettings()
    if find_global_opt:
        scip_settings.MIPGap = 0.0 # must find global optimum
        scip_settings.TimeLimit = t_max_global_opt
    else:
        scip_settings.MIPGap = 1.0 # any feasible solution is acceptable
        scip_settings.TimeLimit = settings.t_max
    scip_settings.VerbLevel = 4 if settings.verbose else 0
    scip_settings.MIPGapAbs = 0.0
    scip_settings.Threads = 1
    scip_settings.FeasibilityTol = settings.eps_prim

    sol = zono.OptSolution()
    Z.optimize_over(P, q, c=c, settings=scip_settings, solution=sol)

    return sol.converged, sol.run_time, sol.z, np.nan

def get_residual(Z, xi, order=2):
    
    # inequalities
    xi_ineq_err = np.zeros(Z.get_nG()) # init
    if Z.is_0_1_form():
        l = np.zeros(Z.get_nG())
    else:
        l = -np.ones(Z.get_nG())
    u = np.ones(Z.get_nG())

    for i in range(Z.get_nGc()):
        if xi[i] < l[i]:
            xi_ineq_err[i] = l[i] - xi[i]
        elif xi[i] > u[i]:
            xi_ineq_err[i] = xi[i] - u[i]
    for i in range(Z.get_nGc(), Z.get_nG()):
        xi_ineq_err[i] = min(np.abs(xi[i] - l[i]), np.abs(u[i] - xi[i]))

    Ax_minus_b = Z.get_A().dot(xi) - Z.get_b()
    return np.linalg.norm(np.hstack( (xi_ineq_err, Ax_minus_b) ), ord=order)

def get_objective(Z, xi, P, q, c):
    x = Z.get_G().dot(xi) + Z.get_c()
    return 0.5*x.dot(P.dot(x)) + q.dot(x) + c


# solver key, legend label, color (Okabe-Ito colorblind-safe palette), marker
SOLVERS = [
    ('admm_fp',  r'ADMM-FP',  '#0072B2', 'x'),  # blue
    ('ofp',      r'OFP',      '#CC79A7', '^'),  # reddish purple
    ('admm',     r'ADMM',     '#009E73', 's'),  # bluish green
    ('gurobi',   r'Gurobi',   '#D55E00', 'o'),  # vermillion
    ('scip',     r'SCIP',     '#56B4E9', 'd'),  # sky blue
]

SOLVE_FNS = {
    'admm_fp': admm_fp_solve,
    'ofp': ofp_solve,
    'admm': admm_solve,
    'gurobi': lambda Z, P, q, c, settings, seed: gurobi_solve(Z, P, q, c, settings, seed, find_global_opt=False),
    'scip': lambda Z, P, q, c, settings, seed: scip_solve(Z, P, q, c, settings, seed, find_global_opt=False),
}


### MAIN ###
n = 100
nGc = 200
nGb = 50
nC = 50
hz_density = 0.1
np.random.seed(0)

settings = zono.OptSettings()
settings.t_max = 1.0
settings.verbose = False

t_arr = {key: [] for key, *_ in SOLVERS}
rel_obj_arr = {key: [] for key, *_ in SOLVERS}
iter_arr = {key: [] for key, *_ in SOLVERS}

n_trials = 100
trial = 0
seed = 0

while trial < n_trials:

    Z = random_hz(n, nGc, nGb, nC, hz_density)

    P = sparse.csc_matrix( (n, n) )
    q = random_vector(n, 1.0)
    c = 0.0

    feas, time_, xi, iters = {}, {}, {}, {}
    for key, *_ in SOLVERS:
        feas[key], time_[key], xi[key], iters[key] = SOLVE_FNS[key](Z, P, q, c, settings, seed=seed)

    global_feas, global_time, global_xi, global_iter = gurobi_solve(Z, P, q, c, settings, seed=seed, find_global_opt=True, t_max_global_opt=10.0)
    global_obj = get_objective(Z, global_xi, P, q, c) if global_feas else np.inf

    obj, rel_obj = {}, {}
    for key, *_ in SOLVERS:
        obj[key] = get_objective(Z, xi[key], P, q, c) if feas[key] else np.inf
        rel_obj[key] = 100.*np.abs(obj[key]-global_obj)/np.abs(obj[key]) if feas[key] else np.inf

    for key, label, *_ in SOLVERS:
        print(f'{label}: feas={feas[key]}, time={time_[key]:.4f}s, computed residual = {get_residual(Z, xi[key]):.4e}, percent suboptimality = {rel_obj[key]:.4e}')

    if any(feas[key] for key, *_ in SOLVERS):
        for key, *_ in SOLVERS:
            if feas[key]:
                t_arr[key].append(time_[key])
                rel_obj_arr[key].append(rel_obj[key])
                iter_arr[key].append(iters[key])
        trial += 1
    seed += 1

# print number of cases where solution is found
print(f'Out of {n_trials} trials:')
for key, label, *_ in SOLVERS:
    print(f'  {label} found feasible solution in {len(t_arr[key])} cases.')

# Save results as JSON
results_dict = {
    'n': n,
    'nGc': nGc,
    'nGb': nGb,
    'nC': nC,
    'hz_density': hz_density,
    'n_trials': n_trials,
    't_arr': t_arr,
    'rel_obj_arr': rel_obj_arr,
    'iter_arr': iter_arr,
}

results_file = 'heuristic_comparison_results.json'
with open(results_file, 'w') as json_file:
    json.dump(results_dict, json_file, indent=4)
print(f'Results written to {results_file}')