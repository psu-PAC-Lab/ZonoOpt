import zonoopt as zono
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import subprocess
import json
import ctypes
import ctypes.util
import multiprocessing as mp
import warnings
from shapely.geometry import Polygon

### flags

# options: 'default', 'admm_fp_multirun', 'heuristic_test'
MODE = 'heuristic_test'

# quick end-to-end check of the heuristic_test pipeline: few seeds, short time limit,
# results written to a separate file so the full dataset isn't overwritten
SMOKE_TEST = True

try:
    import polypartition
except ImportError:
    print('Need to install polypartition from https://github.com/psu-PAC-Lab/polypartition')
    exit()

if MODE == 'heuristic_test' or MODE == 'admm_fp_multirun':
    if MODE == 'heuristic_test':
        try:
            import objective_feasibility_pump as ofp
        except ImportError:
            print('Need to install objective_feasibility_pump from https://github.com/jrobbins11/objective-feasibility-pump')
            exit()
    

### functions
def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

### peak memory measurement
#
# The solvers allocate in C++ and the pybind11 bindings hold the GIL for the duration of the
# solve, so neither tracemalloc nor a Python sampler thread can see the memory a solve uses.
# Instead read the kernel's peak-RSS watermark around the solve: malloc_trim returns heap freed
# during problem construction to the OS, and writing 5 to /proc/self/clear_refs
# (CLEAR_REFS_MMAP_HIWATER_RSS) resets VmHWM so the watermark reflects the solve alone.
#
# Both mechanisms are Linux/glibc specific and neither is required for the rest of this script:
# where they are unavailable a warning is issued, the memory measurement reports nan, and
# everything else runs as before. Where /proc works but malloc_trim doesn't, the measurement
# still works but may under-report a solve that reuses heap freed during problem construction.

def _load_malloc_trim():
    """Return glibc's malloc_trim, or None where it isn't available."""
    try:
        libc = ctypes.CDLL(ctypes.util.find_library('c') or 'libc.so.6')
        return libc.malloc_trim
    except (OSError, AttributeError) as e:
        warnings.warn(f'malloc_trim unavailable ({e}); peak memory measurements may under-report.',
                      RuntimeWarning)
        return None

_malloc_trim = _load_malloc_trim()

def _proc_status_kb(key):
    """Read a memory field (e.g. 'VmRSS:') out of /proc/self/status, in kB."""
    with open('/proc/self/status') as f:
        for line in f:
            if line.startswith(key):
                return int(line.split()[1])
    raise KeyError(key)

def begin_mem_measurement():
    """Reset the peak-RSS watermark. Returns baseline RSS in kB, or None if unsupported."""
    try:
        if _malloc_trim is not None:
            _malloc_trim(0)
        with open('/proc/self/clear_refs', 'w') as f:
            f.write('5\n')
        return _proc_status_kb('VmRSS:')
    except Exception as e:
        warnings.warn(f'peak memory measurement unavailable ({e}); memory results will be nan.',
                      RuntimeWarning)
        return None

def end_mem_measurement(base_kb):
    """Peak memory in MB used since the matching begin_mem_measurement call."""
    if base_kb is None:
        return np.nan
    try:
        return (_proc_status_kb('VmHWM:') - base_kb) / 1024.
    except Exception as e:
        warnings.warn(f'could not read peak memory ({e}); memory result will be nan.', RuntimeWarning)
        return np.nan

# probe once at startup so an unsupported platform warns here rather than mid-run
# (warnings are reported once per call site, so the solve loop below stays quiet)
begin_mem_measurement()

def run_case_in_subprocess(seed, sample_factor, settings, solver_flag):
    """Run one reach-avoid case in a forked child process.

    Isolation matters for the memory measurement: SCIP's block memory pools and Gurobi's
    environment hold on to freed memory, so a second solve by the same solver in the same
    process would reuse resident pages and report far too little. Forking gives every solve a
    clean heap. Solve times are unaffected -- they come from the solvers' own internal timers.

    Returns a dict with converged / run_time / J / iter / mem_mb.
    """

    failed = {'converged': False, 'run_time': np.nan, 'J': np.nan, 'iter': 0, 'mem_mb': np.nan}

    def worker(conn):
        try:
            _, _, _, _, sol, _, mem_mb = reach_avoid_problem(seed, sample_factor, settings, solver_flag)
            conn.send({'converged': sol.converged, 'run_time': sol.run_time,
                       'J': sol.J, 'iter': sol.iter, 'mem_mb': mem_mb})
        except Exception as e:
            conn.send(dict(failed, error=repr(e)))
        finally:
            conn.close()

    # fork so that the (unpicklable) settings object and the already imported zonoopt module
    # are inherited rather than re-created
    ctx = mp.get_context('fork')
    recv_conn, send_conn = ctx.Pipe(False)
    p = ctx.Process(target=worker, args=(send_conn,))
    p.start()
    send_conn.close() # so a dead child gives EOFError instead of blocking forever

    try:
        res = recv_conn.recv()
    except EOFError:
        print(f'{solver_flag} solve died without returning a result (seed {seed}, sample factor {sample_factor})')
        res = failed
    finally:
        recv_conn.close()

    p.join(settings.t_max + 60.)
    if p.is_alive():
        print(f'{solver_flag} solve exceeded its time limit and was terminated (seed {seed}, sample factor {sample_factor})')
        p.terminate()
        p.join()

    if 'error' in res:
        print(f'{solver_flag} solve raised an exception (seed {seed}, sample factor {sample_factor}): {res["error"]}')

    return res

def random_poly(center, n_verts, min_radius, max_radius, rng=None):
    """Makes random polytope"""
    
    # generate random radii
    if rng is not None:
        radii = rng.uniform(min_radius, max_radius, n_verts)
    else:
        radii = np.random.uniform(min_radius, max_radius, n_verts)

    # generate random angles
    if rng is not None:
        angles = rng.uniform(0, 2 * np.pi, n_verts)
    else:
        angles = np.random.uniform(0, 2 * np.pi, n_verts)

    angles = np.sort(angles)

    # generate random vertices
    V = np.zeros((n_verts, 2))
    for i in range(n_verts):
        x = radii[i] * np.cos(angles[i]) + center[0]
        y = radii[i] * np.sin(angles[i]) + center[1]
        V[i] = [x, y]

    # require convexity
    P = Polygon(V)
    P_ch = P.convex_hull

    if P.equals(P_ch): # i.e., is convex
        return P
    else:
        return random_poly(center, n_verts, min_radius, max_radius, rng=rng)
    
def shapely_poly_2_vrep(P):
    """Convert a Shapely Polygon to a list of vertices.
    P : shapely Polygon"""
    
    # get vertices
    V = np.array(P.exterior.coords.xy).transpose()
    V = V[0:-1,:] # remove duplicate vertex

    return V

def zono_planning_prob(x0, xr, A, B, Q, R, QN, N, X, U, XN, X_map, solver_flag, settings=zono.OptSettings(), sol=zono.OptSolution()):
    """Build and solve MPC problem using zonoopt set operations.
    x0 = initial state
    xr = reference trajectory (N+1 x n numpy array)
    A = state transition matrix
    B = control matrix
    Q = state cost matrix
    R = control cost matrix
    QN = terminal state cost matrix
    N = prediction horizon
    X = state constraints
    U = control constraint set
    XN = terminal state constraint
    X_map = map constraints (obs avoidance)
    settings = zonoopt.OptSettings object
    sol = zonoopt.OptSolution object
    """

    # dims
    nx = A.shape[1]
    nu = B.shape[1]

    # index tracking
    idx = 0
    idx_x = []
    idx_u = []

    # initial state
    Z = zono.Point(x0)
    idx_x.append([j for j in range(idx, idx+nx)])
    idx += nx

    # init cost
    P = Q
    q = np.zeros(nx)
    c = 0.0 

    for k in range(N):

        # control
        Z = zono.cartesian_product(Z, U)

        idx_u.append([j for j in range(idx, idx+nu)])
        idx += nu
        
        # dynamics
        nZ = Z.get_n()
        Z = zono.cartesian_product(Z, X)
        ABmI = sp.hstack((sp.csc_matrix((nx, nZ-(nx+nu))), A, B, -sp.eye(nx)))
        Z = zono.intersection(Z, zono.Point(np.zeros(nx)), ABmI)

        idx_x.append([j for j in range(idx, idx+nx)])
        idx += nx

        # map
        Z = zono.intersection_over_dims(Z, X_map, idx_x[-1][0:2])

        # terminal state constraint
        if k == N-1:
            Z = zono.intersection_over_dims(Z, XN, idx_x[-1])

        # cost
        if k == N-1:
            P = sp.block_diag((P, R, QN))
            q = np.hstack([q, np.zeros(nu), -QN.dot(xr[k])])
            c += 0.5*xr[k].dot(QN.dot(xr[k]))
        else:
            P = sp.block_diag((P, R, Q))
            q = np.hstack([q, np.zeros(nu), -Q.dot(xr[k])])
            c += 0.5*xr[k].dot(Q.dot(xr[k]))

    # solve
    # memory measurement brackets the solve call only, so it excludes the set operations above
    if solver_flag == 'zonoopt':
        mem_base = begin_mem_measurement()
        xopt = Z.optimize_over(P, q, c=c, settings=settings, solution=sol)
        mem_mb = end_mem_measurement(mem_base)

    elif solver_flag == 'gurobi':

        # solve using ZonoOpt's built-in Gurobi solver support
        gurobi_settings = zono.GurobiSettings()
        gurobi_settings.MIPGap = settings.eps_r
        gurobi_settings.MIPGapAbs = settings.eps_a
        gurobi_settings.Threads = settings.n_threads_bnb + settings.n_threads_admm_fp + 1  # +1 for the main thread
        gurobi_settings.OutputFlag = 1 if settings.verbose else 0
        gurobi_settings.TimeLimit = settings.t_max
        gurobi_settings.FeasibilityTol = settings.eps_prim

        mem_base = begin_mem_measurement()
        xopt = Z.optimize_over(P, q, c=c, settings=gurobi_settings, solution=sol)
        mem_mb = end_mem_measurement(mem_base)

        if not sol.converged:
            print('Gurobi did not find a feasible solution')
            sol.infeasible = True
            return None, None, None, mem_mb

    elif solver_flag == 'scip':

        # solve using ZonoOpt's built-in SCIP solver support
        scip_settings = zono.SCIPSettings()
        scip_settings.MIPGap = settings.eps_r
        scip_settings.MIPGapAbs = settings.eps_a
        scip_settings.Threads = settings.n_threads_bnb + settings.n_threads_admm_fp + 1  # +1 for the main thread
        scip_settings.VerbLevel = 4 if settings.verbose else 0
        scip_settings.TimeLimit = settings.t_max
        scip_settings.FeasibilityTol = settings.eps_prim

        mem_base = begin_mem_measurement()
        xopt = Z.optimize_over(P, q, c=c, settings=scip_settings, solution=sol)
        mem_mb = end_mem_measurement(mem_base)

        if not sol.converged:
            print('SCIP did not find a feasible solution')
            sol.infeasible = True
            return None, None, None, mem_mb

    elif solver_flag == 'ofp':
        # settings
        ofp_settings = ofp.OFP_Settings()
        ofp_settings.t_max = settings.t_max
        ofp_settings.max_iter = 1000000 # very large number
        ofp_settings.buffer_size = 20
        ofp_settings.verbose = settings.verbose
        ofp_settings.verbosity_interval = settings.verbosity_interval
        ofp_settings.max_restarts = 100000 # very large number
        ofp_settings.rng_seed = settings.rng_seed
        ofp_settings.inf_norm_conv = True
        ofp_settings.tol = settings.eps_prim

        # problem data
        q_tilde = Z.get_G().transpose().dot(P*Z.get_c() + q)
        bins = [i for i in range(Z.get_nG() - Z.get_nGb(), Z.get_nG())]

        # solve
        if not Z.is_0_1_form():
            Z.convert_form()
        solver = ofp.OFP_Solver(q_tilde, Z.get_A(), Z.get_b(), Z.get_b(), np.zeros(Z.get_nG()), np.ones(Z.get_nG()), bins, settings=ofp_settings)
        mem_base = begin_mem_measurement()
        success = solver.solve()
        mem_mb = end_mem_measurement(mem_base)

        if not success:
            print('OFP did not find a feasible solution')
            info = solver.get_info()
            print(f'iter: {info.iter}')
            print(f'restarts: {info.restarts}')
            print(f'perturbations: {info.perturbations}')

            # get residual
            res = 0.0
            Ax_m_b = Z.get_A().dot(solver.get_solution()) - Z.get_b()
            for vi in Ax_m_b:
                res += np.abs(vi)
            for xi_i in solver.get_solution():
                if xi_i < 0.0:
                    res -= xi_i
                elif xi_i > 1.0:
                    res += xi_i - 1.0
            print(f'constraint residual: {res}')

        # logging
        info = solver.get_info()
        sol.run_time = info.runtime
        sol.iter = info.iter
        sol.startup_time = 0.
        sol.infeasible = False
        sol.J = info.objective
        sol.converged = info.feasible

        # solution
        xi_opt = solver.get_solution()
        xopt = Z.get_G()*xi_opt + Z.get_c()

    else:
        raise ValueError(f'Unknown solver_flag: {solver_flag}')


    # state, input trajectory
    x_traj = []
    u_traj = []
    for idx in idx_x:
        x_traj.append(xopt[idx])
    for idx in idx_u:
        u_traj.append(xopt[idx])

    # return
    return x_traj, u_traj, Z, mem_mb


### make map
def reach_avoid_problem(seed, sample_factor, settings, solver_flag):
    
    # map parameters
    min_radius = 1.0
    max_radius = 2.0
    x_min = 0.
    x_max = 10.
    y_min = -5.
    y_max = 5.
    n_obs = 3
    n_sides = 5
    dw = (max_radius+min_radius)/2.

    # boundary
    Vbdy = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])
    Pbdy = Polygon(Vbdy)

    # random obstacles
    rng = np.random.default_rng(seed)
    Pobs = []
    for i in range(n_obs):

        intersects_bdy = True
        intersects_polys = True
        while intersects_bdy or intersects_polys:
        
            cx = rng.uniform(x_min+dw, x_max-dw)
            cy = rng.uniform(y_min+dw, y_max-dw)
            poly = random_poly(np.array([cx, cy]), n_sides, min_radius, max_radius, rng=rng)

            intersects_bdy = poly.bounds[0] < x_min or poly.bounds[1] < y_min or poly.bounds[2] > x_max or poly.bounds[3] > y_max
            intersects_polys = False
            for P in Pobs:
                if P.intersects(poly):
                    intersects_polys = True
                    break

        Pobs.append(poly)

    Vobs = [shapely_poly_2_vrep(P) for P in Pobs]

    # partition
    Vpart = polypartition.hertel_mehlhorn(Vbdy, Vobs)

    # convert to hybzono
    Z_map = zono.vrep_2_hybzono([Vpoly for Vpoly in Vpart])

    ### motion planning problem ###

    # time step
    dt = 2.0/sample_factor

    # 2D double integrator dynamics
    # x = [x, y, xdot, y_dot]
    # u = [x_ddot, y_ddot]
    A = sp.csc_matrix(np.array(
                [[1., 0., dt, 0.],
                [0., 1., 0., dt],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]))
    B = sp.csc_matrix(np.array(
                [[0.5*dt**2, 0.],
                [0., 0.5*dt**2],
                [dt, 0.],
                [0., dt]]))
    nx = 4
    nu = 2

    # state feasible set

    # velocity constraints
    v_max = 1.0
    Xv = zono.make_regular_zono_2D(v_max, 4)

    # full state constraints
    X_map_outer = zono.interval_2_zono(zono.Box(np.array([x_min, y_min]), np.array([x_max, y_max])))
    X = zono.cartesian_product(X_map_outer, Xv)
    # X = zono.cartesian_product(Z_map, Xv)

    # input feasible set

    # turn rate constraints
    om_max = 90*np.pi/180
    v_min = 0.1 # fictitious, om_max not enforced below this
    U = zono.make_regular_zono_2D(om_max*v_min, 4)

    # MPC horizon
    N = 10*sample_factor

    # cost function matrices
    # Q = sp.csc_matrix((nx, nx))
    QN = sp.diags([1.0, 1.0, 0.0, 0.0])
    Q = 0.1*(1./N)*QN
    R = (10/N)*sp.eye(2)

    # initial condition
    x = np.array([x_min+0.1, (y_min+y_max)/2., 0.1, 0.])

    # reference
    x_target = np.array([x_max, (y_min+y_max)/2., 0., 0.])
    xr = []
    for i in range(N):
        xr.append(x_target)

    # terminal constraint
    dXN_p = zono.make_regular_zono_2D(1.0, 6)
    dXn_v = zono.make_regular_zono_2D(0.01, 6)
    dXN = zono.cartesian_product(dXN_p, dXn_v)  
    XN = zono.minkowski_sum(zono.Point(x_target), dXN)

    # solve motion planning problem
    sol = zono.OptSolution()
    x_traj, u_traj, Z_mpc, mem_mb = zono_planning_prob(x, xr, A, B, Q, R, QN, N, X, U, XN, Z_map, solver_flag, settings=settings, sol=sol)
    return [x_traj, u_traj, Z_mpc, Z_map, sol, XN, mem_mb]

### main execution
if MODE == 'admm_fp_multirun':

    # settings
    settings = zono.OptSettings()
    settings.t_max = 10.0
    settings.verbose = False
    settings.verbosity_interval = 1000

    # parameters to do single-threaded, one-shot ADMM-FP solution
    settings.single_threaded_admm_fp = True
    settings.polish = False
    settings.enable_rng_seed = True
    settings.k_max_admm_fp_ph2 = int(1e7) # massive number, just let max time be limiting factor

    # build and solve motion planning problem
    seed = 13
    sample_factor = 4

    # number of tests to run
    n_tests = 25

    # run tests
    sol_time_arr = []
    objective_arr = []
    x_traj_arr = []
    for rng_seed in range(n_tests):
        print(f'Running test {rng_seed+1} / {n_tests}')

        settings.rng_seed = rng_seed
        x_traj, u_traj, Z_mpc, Z_map, sol, XN, _ = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')
        if sol.converged:
            x_traj_arr.append(x_traj)
            sol_time_arr.append(sol.run_time)
            objective_arr.append(sol.J)

    # get global optimum for comparison with Gurobi
    settings.n_threads_bnb = 1
    settings.n_threads_admm_fp = 0
    settings.eps_r = 0.0
    settings.eps_a = 0.0
    x_traj_global, u_traj, Z_mpc, Z_map, sol, XN, _ = reach_avoid_problem(seed, sample_factor, settings, 'gurobi')
    objective_global = sol.J

    J_rel = lambda J : np.abs(J - objective_global) / np.abs(J)

    # print results
    print(f'median solution time = {np.median(sol_time_arr)} s, min sol time = {np.min(sol_time_arr)} s, max sol time = {np.max(sol_time_arr)} s')
    print(f'median sub-optimality = {J_rel(np.median(objective_arr))*100}%, min sub-optimality = {J_rel(np.min(objective_arr))*100}%, max sub-optimality = {J_rel(np.max(objective_arr))*100}%')
    print(f'number of cases converged = {len(sol_time_arr)} / {n_tests}')
    
    ### plot ###
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
    figsize = (245.71 * inches_per_pt, 1*245.71 * inches_per_pt)  # Convert pt to inches

    # plot
    with plt.rc_context(rc_context):

        # figure
        fig = plt.figure(constrained_layout=True, figsize=figsize)
        ax = fig.add_subplot(111)

        # plot map and terminal constraint
        zono.plot(Z_map, ax=ax, color='b', alpha=0.25, edgecolor='b', linewidth=1.0)
        zono.plot(zono.intersection(zono.project_onto_dims(XN, [0,1]), Z_map), ax=ax, color='g', alpha=0.5, edgecolor='g', linewidth=1.0)

        # trajectory solution
        for x_traj in x_traj_arr:
            x_vec = np.array(x_traj)
            ax.plot(x_vec[:,0], x_vec[:,1], '.r', markersize=2.5, label=r'ADMM-FP')
        
        # global optimum
        x_vec_global = np.array(x_traj_global)
        ax.plot(x_vec_global[:,0], x_vec_global[:,1], '.k', markersize=2.5, label=r'Global Optimum')

        # axes
        ax.set_xlabel(r'$x$ [m]')
        ax.set_ylabel(r'$y$ [m]')
        ax.axis('equal')
        ax.grid(alpha=0.2)

        # make legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=textwidth_pt)

        # save
        if is_latex_installed():
            plt.savefig('reach_avoid_traj.pgf')

        plt.show()

elif MODE == 'heuristic_test':

    # settings
    settings = zono.OptSettings()
    settings.t_max = 30.
    settings.verbose = False

    # parameters to do single-threaded, one-shot ADMM-FP solution
    settings.single_threaded_admm_fp = True
    settings.polish = False
    settings.use_interval_contractor = False
    settings.rng_seed = 0
    settings.enable_rng_seed = True

    # settings needed for Gurobi
    settings.eps_r = 1.0
    settings.eps_a = 0.0
    settings.n_threads_bnb = 0
    settings.n_threads_admm_fp = 0
    
    # loop through seeds, trials, sample factors
    n_seeds = 100
    sample_factor_arr = [2*i-1 for i in range(1,6)]
    results_file = 'reach_avoid_results.json'

    # abbreviated run for checking the pipeline end to end
    if SMOKE_TEST:
        n_seeds = 5
        settings.t_max = 5.
        results_file = 'reach_avoid_results_smoke.json'
        print(f'SMOKE_TEST enabled: {n_seeds} seeds, sample factors {sample_factor_arr}, '
              f't_max = {settings.t_max} s, writing to {results_file}')

    # configurations for different ADMM tests
    settings_admm_fp = settings.copy()
    settings_admm_fp.enable_perturb_admm_fp = True
    settings_admm_fp.enable_restart_admm_fp = True
    settings_admm_fp.k_max_admm_fp_ph1 = 10000 # default
    settings_admm_fp.k_max_admm_fp_ph2 = int(1e7) # large number

    settings_admm = settings.copy()
    settings_admm.enable_perturb_admm_fp = False
    settings_admm.enable_restart_admm_fp = False
    settings_admm.k_max_admm_fp_ph1 = int(1e7) # large number
    settings_admm.k_max_admm_fp_ph2 = 0

    settings_admm_rnd = settings_admm.copy()
    settings_admm_rnd.k_max_admm_fp_ph1 = 1000 # small number
    settings_admm_rnd.polish = True

    n_feas_admm_fp_arr = []
    n_feas_admm_arr = []
    n_feas_admm_rnd_arr = []
    n_feas_gurobi_arr = []
    n_feas_scip_arr = []
    n_feas_ofp_arr = []
    t_admm_fp_arr = []
    t_admm_arr = []
    t_admm_rnd_arr = []
    t_gurobi_arr = []
    t_scip_arr = []
    t_ofp_arr = []
    m_admm_fp_arr = []
    m_admm_arr = []
    m_admm_rnd_arr = []
    m_gurobi_arr = []
    m_scip_arr = []
    m_ofp_arr = []

    for sample_factor in sample_factor_arr:

        # reset counters
        n_feas_admm = 0
        n_feas_admm_fp = 0
        n_feas_admm_rnd = 0
        n_feas_gurobi = 0
        n_feas_scip = 0
        n_feas_ofp = 0
        t_admm_fp = []
        t_admm = []
        t_admm_rnd = []
        t_gurobi = []
        t_scip = []
        t_ofp = []
        m_admm_fp = []
        m_admm = []
        m_admm_rnd = []
        m_gurobi = []
        m_scip = []
        m_ofp = []

        for seed in range(n_seeds):
            
            print(f'Running sample factor {sample_factor}, seed {seed}')

            # ADMM-FP
            res = run_case_in_subprocess(seed, sample_factor, settings_admm_fp, 'zonoopt')

            if res['converged']:
                n_feas_admm_fp += 1
                t_admm_fp.append(res['run_time'])
                m_admm_fp.append(res['mem_mb'])
            else:
                t_admm_fp.append(np.nan)
                m_admm_fp.append(np.nan)

            # ADMM
            res = run_case_in_subprocess(seed, sample_factor, settings_admm, 'zonoopt')

            if res['converged']:
                n_feas_admm += 1
                t_admm.append(res['run_time'])
                m_admm.append(res['mem_mb'])
            else:
                t_admm.append(np.nan)
                m_admm.append(np.nan)

            # ADMM with reduced iterations and rounding
            res = run_case_in_subprocess(seed, sample_factor, settings_admm_rnd, 'zonoopt')
            
            if res['converged']:
                n_feas_admm_rnd += 1
                t_admm_rnd.append(res['run_time'])
                m_admm_rnd.append(res['mem_mb'])
            else:
                t_admm_rnd.append(np.nan)
                m_admm_rnd.append(np.nan)

            # run with Gurobi for comparison
            res = run_case_in_subprocess(seed, sample_factor, settings, 'gurobi')

            if res['converged']:
                n_feas_gurobi += 1
                t_gurobi.append(res['run_time'])
                m_gurobi.append(res['mem_mb'])
            else:
                t_gurobi.append(np.nan)
                m_gurobi.append(np.nan)

            # run with SCIP for comparison
            res = run_case_in_subprocess(seed, sample_factor, settings, 'scip')

            if res['converged']:
                n_feas_scip += 1
                t_scip.append(res['run_time'])
                m_scip.append(res['mem_mb'])
            else:
                t_scip.append(np.nan)
                m_scip.append(np.nan)

            # run with OFP for comparison
            res = run_case_in_subprocess(seed, sample_factor, settings, 'ofp')

            if res['converged']:
                n_feas_ofp += 1
                t_ofp.append(res['run_time'])
                m_ofp.append(res['mem_mb'])
            else:
                t_ofp.append(np.nan)
                m_ofp.append(np.nan)


        # log results
        n_feas_admm_fp_arr.append(n_feas_admm_fp)
        n_feas_admm_arr.append(n_feas_admm)
        n_feas_admm_rnd_arr.append(n_feas_admm_rnd)
        n_feas_gurobi_arr.append(n_feas_gurobi)
        n_feas_scip_arr.append(n_feas_scip)
        n_feas_ofp_arr.append(n_feas_ofp)
        t_admm_fp_arr.append(t_admm_fp)
        t_admm_arr.append(t_admm)
        t_admm_rnd_arr.append(t_admm_rnd)
        t_gurobi_arr.append(t_gurobi)
        t_scip_arr.append(t_scip)
        t_ofp_arr.append(t_ofp)
        m_admm_fp_arr.append(m_admm_fp)
        m_admm_arr.append(m_admm)
        m_admm_rnd_arr.append(m_admm_rnd)
        m_gurobi_arr.append(m_gurobi)
        m_scip_arr.append(m_scip)
        m_ofp_arr.append(m_ofp)

        # print results
        # nanmean over an all-nan sample factor warns rather than failing, so silence that case
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            for label, n_feas, t, m in [('ADMM-FP', n_feas_admm_fp, t_admm_fp, m_admm_fp),
                                        ('ADMM', n_feas_admm, t_admm, m_admm),
                                        ('ADMM-RND', n_feas_admm_rnd, t_admm_rnd, m_admm_rnd),
                                        ('Gurobi', n_feas_gurobi, t_gurobi, m_gurobi),
                                        ('SCIP', n_feas_scip, t_scip, m_scip),
                                        ('OFP', n_feas_ofp, t_ofp, m_ofp)]:
                print(f'{label} feasible in {n_feas} / {n_seeds} cases -> {n_feas/(n_seeds)*100:.2f}%, '
                      f'average solution time = {np.nanmean(t) if len(t) > 0 else np.nan} s, '
                      f'median peak memory = {np.nanmedian(m) if len(m) > 0 else np.nan} MB')

    # Save results as JSON
    results_dict = {
        'sample_factor_arr': sample_factor_arr,
        'n_seeds': n_seeds,
        'n_feas_admm_fp_arr': n_feas_admm_fp_arr,
        'n_feas_admm_arr': n_feas_admm_arr,
        'n_feas_admm_rnd_arr': n_feas_admm_rnd_arr,
        'n_feas_gurobi_arr': n_feas_gurobi_arr,
        'n_feas_scip_arr': n_feas_scip_arr,
        'n_feas_ofp_arr': n_feas_ofp_arr,
        't_admm_fp_arr': t_admm_fp_arr,
        't_admm_arr': t_admm_arr,
        't_admm_rnd_arr': t_admm_rnd_arr,
        't_gurobi_arr': t_gurobi_arr,
        't_scip_arr': t_scip_arr,
        't_ofp_arr': t_ofp_arr,
        'm_admm_fp_arr': m_admm_fp_arr,
        'm_admm_arr': m_admm_arr,
        'm_admm_rnd_arr': m_admm_rnd_arr,
        'm_gurobi_arr': m_gurobi_arr,
        'm_scip_arr': m_scip_arr,
        'm_ofp_arr': m_ofp_arr
    }

    with open(results_file, 'w') as json_file:
        json.dump(results_dict, json_file, indent=4)
    print(f'Results written to {results_file}')


else:

    # settings
    settings = zono.OptSettings()
    settings.eps_r = 1.0
    settings.rho = 10.
    settings.t_max = 30.0
    settings.verbose = True
    settings.verbosity_interval = 1000
    settings.single_threaded_admm_fp = True
    
    # build and solve motion planning problem
    seed = 13
    sample_factor = 4
    x_traj, u_traj, Z_mpc, Z_map, sol, XN, _ = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')

    print(f'solution time = {sol.run_time} s, startup time = {sol.startup_time} s, number of nodes evaluated = {sol.iter}')
    print(f'Z_mpc number of continuous variables = {Z_mpc.get_nGc()}, number of binary variables = {Z_mpc.get_nGb()}, number of constraints = {Z_mpc.get_nC()}')

    ### plot ###
    textwidth_pt = 12.
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
    figsize = (245.71 * inches_per_pt, 1*245.71 * inches_per_pt)  # Convert pt to inches

    # plot
    with plt.rc_context(rc_context):

        # figure
        fig = plt.figure(constrained_layout=True, figsize=figsize)
        ax = fig.add_subplot(111)

        # time-varying state constraints (spatio-temporal corridor)
        zono.plot(Z_map, ax=ax, color='b', alpha=0.25, edgecolor='b', linewidth=1.0)
        zono.plot(zono.intersection(zono.project_onto_dims(XN, [0,1]), Z_map), ax=ax, color='g', alpha=0.5, edgecolor='g', linewidth=1.0)

        # solution
        x_vec = np.array(x_traj)
        ax.plot(x_vec[:,0], x_vec[:,1], '.r', markersize=5)

        # axes
        ax.set_xlabel('$x$ [m]')
        ax.set_ylabel('$y$ [m]')
        ax.axis('equal')
        ax.grid(alpha=0.2)

        plt.show()