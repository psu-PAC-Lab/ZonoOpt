import numpy as np
from scipy import sparse
import zonoopt as zono
try:
    import objective_feasibility_pump as ofp
except ImportError as e:
    print("Need to pip install OFP from https://github.com/jrobbins11/objective-feasibility-pump")
    raise e
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import subprocess

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

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

t_admm_fp_arr = []
t_admm_arr = []
t_ofp_arr = []
t_gurobi_arr = []
t_scip_arr = []
rel_obj_admm_fp_arr = []
rel_obj_admm_arr = []
rel_obj_ofp_arr = []
rel_obj_gurobi_arr = []
rel_obj_scip_arr = []
iter_admm_fp_arr = []
iter_admm_arr = []
iter_ofp_arr = []
iter_gurobi_arr = []
iter_scip_arr = []

n_trials = 100
trial = 0
seed = 0

while trial < n_trials:

    Z = random_hz(n, nGc, nGb, nC, hz_density)

    P = sparse.csc_matrix( (n, n) )
    q = random_vector(n, 1.0)
    c = 0.0

    admm_fp_feas, admm_fp_time, admm_fp_xi, admm_fp_iter = admm_fp_solve(Z, P, q, c, settings, seed=seed)
    admm_feas, admm_time, admm_xi, admm_iter = admm_solve(Z, P, q, c, settings, seed=seed)
    ofp_feas, ofp_time, ofp_xi, ofp_iter = ofp_solve(Z, P, q, c, settings, seed=seed)
    gurobi_feas, gurobi_time, gurobi_xi, gurobi_iter = gurobi_solve(Z, P, q, c, settings, seed=seed, find_global_opt=False)
    scip_feas, scip_time, scip_xi, scip_iter = scip_solve(Z, P, q, c, settings, seed=seed, find_global_opt=False)
    global_feas, global_time, global_xi, global_iter = gurobi_solve(Z, P, q, c, settings, seed=seed, find_global_opt=True, t_max_global_opt=10.0)

    admm_fp_obj = get_objective(Z, admm_fp_xi, P, q, c) if admm_fp_feas else np.inf
    admm_obj = get_objective(Z, admm_xi, P, q, c) if admm_feas else np.inf
    ofp_obj = get_objective(Z, ofp_xi, P, q, c) if ofp_feas else np.inf
    gurobi_obj = get_objective(Z, gurobi_xi, P, q, c) if gurobi_feas else np.inf
    scip_obj = get_objective(Z, scip_xi, P, q, c) if scip_feas else np.inf
    global_obj = get_objective(Z, global_xi, P, q, c) if global_feas else np.inf

    admm_fp_rel_obj = 100.*np.abs(admm_fp_obj-global_obj)/np.abs(admm_fp_obj) if admm_fp_feas else np.inf
    admm_rel_obj = 100.*np.abs(admm_obj-global_obj)/np.abs(admm_obj) if admm_feas else np.inf
    ofp_rel_obj = 100.*np.abs(ofp_obj-global_obj)/np.abs(ofp_obj) if ofp_feas else np.inf
    gurobi_rel_obj = 100.*np.abs(gurobi_obj-global_obj)/np.abs(gurobi_obj) if gurobi_feas else np.inf
    scip_rel_obj = 100.*np.abs(scip_obj-global_obj)/np.abs(scip_obj) if scip_feas else np.inf

    print(f'ADMM-FP:   feas={admm_fp_feas}, time={admm_fp_time:.4f}s, computed residual = {get_residual(Z, admm_fp_xi):.4e}, percent suboptimality = {admm_fp_rel_obj:.4e}')
    print(f'ADMM:      feas={admm_feas}, time={admm_time:.4f}s, computed residual = {get_residual(Z, admm_xi):.4e}, percent suboptimality = {admm_rel_obj:.4e}')
    print(f'OFP:      feas={ofp_feas}, time={ofp_time:.4f}s, computed residual = {get_residual(Z, ofp_xi):.4e}, percent suboptimality = {ofp_rel_obj:.4e}')
    print(f'Gurobi:   feas={gurobi_feas}, time={gurobi_time:.4f}s, computed residual = {get_residual(Z, gurobi_xi):.4e}, percent suboptimality = {gurobi_rel_obj:.4e}')
    print(f'SCIP:     feas={scip_feas}, time={scip_time:.4f}s, computed residual = {get_residual(Z, scip_xi):.4e}, percent suboptimality = {scip_rel_obj:.4e}')

    if admm_fp_feas or admm_feas or ofp_feas or gurobi_feas or scip_feas:
        if admm_fp_feas:
            t_admm_fp_arr.append(admm_fp_time)
            rel_obj_admm_fp_arr.append(admm_fp_rel_obj)
            iter_admm_fp_arr.append(admm_fp_iter)
        if admm_feas:
            t_admm_arr.append(admm_time)
            rel_obj_admm_arr.append(admm_rel_obj)
            iter_admm_arr.append(admm_iter)
        if ofp_feas:
            t_ofp_arr.append(ofp_time)
            rel_obj_ofp_arr.append(ofp_rel_obj)
            iter_ofp_arr.append(ofp_iter)
        if gurobi_feas:
            t_gurobi_arr.append(gurobi_time)
            rel_obj_gurobi_arr.append(gurobi_rel_obj)
            iter_gurobi_arr.append(gurobi_iter)
        if scip_feas:
            t_scip_arr.append(scip_time)
            rel_obj_scip_arr.append(scip_rel_obj)
            iter_scip_arr.append(scip_iter)
        trial += 1
    seed += 1

# print number of cases where solution is found
print(f'Out of {n_trials} trials:')
print(f'  ADMM-FP found feasible solution in {len(t_admm_fp_arr)} cases.')
print(f'  ADMM found feasible solution in {len(t_admm_arr)} cases.')
print(f'  OFP found feasible solution in {len(t_ofp_arr)} cases.')
print(f'  Gurobi found feasible solution in {len(t_gurobi_arr)} cases.')
print(f'  SCIP found feasible solution in {len(t_scip_arr)} cases.')

# plot results
textwidth_pt = 10.
if is_latex_installed():
    rc_context = {
        "text.usetex": True,
        "font.size": textwidth_pt,
        "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
    }
else:
    print("LaTeX not installed, using default font.")
    rc_context = {
        "font.size": textwidth_pt,
    }


inches_per_pt = 1 / 72.27

with plt.rc_context(rc_context):

    data_labels = [r'ADMM-FP', r'OFP', r'ADMM', r'Gurobi', r'SCIP']
    colors = ['b', 'm', 'g', 'r', 'c']

    # solution time / rate of solution comparison plot
    figwidth_pt = 245.71
    figsize = (figwidth_pt * inches_per_pt, 1.3*figwidth_pt * inches_per_pt)  # Convert pt to inches
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = fig.add_gridspec(ncols=1, nrows=4, height_ratios=[1, 2, 2, 2])

    ax = fig.add_subplot(gs[0,0])
    ax.bar(range(5),
           [len(t_admm_fp_arr)*(100./n_trials), len(t_ofp_arr)*(100./n_trials), len(t_admm_arr)*(100./n_trials), len(t_gurobi_arr)*(100./n_trials), len(t_scip_arr)*(100./n_trials)],
           color=colors,
           alpha=0.5)
    ax.set_xticklabels([])
    ax.set_ylabel(r'[\%]', fontsize=textwidth_pt)
    ax.set_title(r'Percent trials with feasible solution found', fontsize=textwidth_pt)
    ax.grid(axis='y', which='major', alpha=0.2)


    # suboptimality of found solution
    ax = fig.add_subplot(gs[1,0])

    all_data = [rel_obj_admm_fp_arr, rel_obj_ofp_arr, rel_obj_admm_arr, rel_obj_gurobi_arr, rel_obj_scip_arr]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.) # whiskers cover all data
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Suboptimality of found solution', fontsize=textwidth_pt)
    ax.set_ylabel(r'[\%]', fontsize=textwidth_pt)

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.set_xticklabels([])


    # time to find a feasible solution
    ax = fig.add_subplot(gs[2,0])

    all_data = [t_admm_fp_arr, t_ofp_arr, t_admm_arr, t_gurobi_arr, t_scip_arr]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.) # whiskers cover all data
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Time to find a feasible solution', fontsize=textwidth_pt)
    ax.set_ylabel(r'[sec]', fontsize=textwidth_pt)
    ax.set_yscale('log')

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.set_xticklabels([])


    # number of iterations to find a feasible solution
    ax = fig.add_subplot(gs[3,0])

    all_data = [iter_admm_fp_arr, iter_ofp_arr, iter_admm_arr, iter_gurobi_arr, iter_scip_arr]
    bp = ax.boxplot(
        all_data,
        patch_artist=True,
        medianprops={'color': 'black'},
        showfliers=True,
        sym='k.',
        whis=(0., 100.), # whiskers cover all data
        tick_labels=data_labels
    )

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        box.set_alpha(0.5)

    ax.set_title(r'Number of iterations', fontsize=textwidth_pt)
    ax.set_ylabel(r'[num]', fontsize=textwidth_pt)
    ax.set_yscale('log')

    ax.grid(axis='y', which='major', alpha=0.2)

    ax.tick_params(axis='x', which='minor', bottom=False)

    if is_latex_installed():
        plt.savefig('heuristic_comparison.pgf')

    plt.show()