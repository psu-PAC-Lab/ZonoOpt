import zonoopt as zono
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import subprocess
import json
from shapely.geometry import Polygon

### flags

# options: 'default', 'admm_fp_multirun', 'heuristic_test'
MODE = 'default'

try:
    import polypartition
except ImportError:
    print('Need to install polypartition from https://github.com/psu-PAC-Lab/polypartition')
    exit()

if MODE == 'heuristic_test' or MODE == 'admm_fp_multirun':
    import gurobipy as gp
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
    if solver_flag == 'zonoopt':
        xopt = Z.optimize_over(P, q, c=c, settings=settings, solution=sol)
    
    elif solver_flag == 'gurobi':

        # create model
        prob = gp.Model()

        # settings
        prob.Params.MIPGap = settings.eps_r
        prob.Params.MIPGapAbs = settings.eps_a
        prob.Params.Threads = settings.n_threads_bnb + settings.n_threads_admm_fp + 1  # +1 for the main thread
        prob.Params.OutputFlag = 1 if settings.verbose else 0
        prob.Params.TimeLimit = settings.t_max
        prob.Params.FeasibilityTol = settings.eps_prim

        # need 0-1 form
        if not Z.is_0_1_form():
            Z.convert_form()

        # build gurobi model
        P_tilde = Z.get_G().transpose().dot(P.dot(Z.get_G()))
        q_tilde = Z.get_G().transpose().dot(P*Z.get_c() + q)
        c_tilde = 0.5*Z.get_c().dot(P.dot(Z.get_c())) + q.dot(Z.get_c()) + c

        # add variables
        x_gurobi = prob.addMVar(Z.get_nG(), lb=np.zeros(Z.get_nG()), ub=np.ones(Z.get_nG()))
        xc = x_gurobi[0:Z.get_nGc()]
        xb = x_gurobi[Z.get_nGc():Z.get_nG()]
        xc.vtype = gp.GRB.CONTINUOUS
        xb.vtype = gp.GRB.BINARY

        # add constraints
        prob.addConstr(Z.get_Ac().dot(xc) + Z.get_Ab().dot(xb) == Z.get_b())

        # add objective
        prob.setMObjective(0.5*P_tilde, q_tilde, c_tilde, sense=gp.GRB.MINIMIZE)

        # optimize
        prob.optimize()

        if prob.Status != 2:
            print('Gurobi did not find a feasible solution')
            sol.infeasible = True
            return None, None, None

        # logging
        sol.run_time = prob.Runtime
        sol.iter = prob.BarIterCount
        sol.startup_time = 0.0
        sol.infeasible = prob.Status != 2
        sol.J = prob.ObjVal
        sol.converged = prob.Status == 2

        # solution
        xi_opt = np.array([x.X for x in prob.getVars()])
        xopt = Z.get_G()*xi_opt + Z.get_c()

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
        success = solver.solve()

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
    return x_traj, u_traj, Z


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
    x_traj, u_traj, Z_mpc = zono_planning_prob(x, xr, A, B, Q, R, QN, N, X, U, XN, Z_map, solver_flag, settings=settings, sol=sol)
    return [x_traj, u_traj, Z_mpc, Z_map, sol, XN]

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
        x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')
        if sol.converged:
            x_traj_arr.append(x_traj)
            sol_time_arr.append(sol.run_time)
            objective_arr.append(sol.J)

    # get global optimum for comparison with Gurobi
    settings.n_threads_bnb = 1
    settings.n_threads_admm_fp = 0
    settings.eps_r = 0.0
    settings.eps_a = 0.0
    x_traj_global, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'gurobi')
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
    sample_factor_arr = [i for i in range(1,11)]
    
    n_feas_admm_fp_arr = []
    n_feas_admm_arr = []
    n_feas_gurobi_arr = []
    n_feas_ofp_arr = []
    t_admm_fp_arr = []
    t_admm_arr = []
    t_gurobi_arr = []
    t_ofp_arr = []

    for sample_factor in sample_factor_arr:

        # reset counters
        n_feas_admm = 0
        n_feas_admm_fp = 0
        n_feas_gurobi = 0
        n_feas_ofp = 0
        t_admm_fp = []
        t_admm = []
        t_gurobi = []
        t_ofp = []

        for seed in range(n_seeds):
            
            print(f'Running sample factor {sample_factor}, seed {seed}')

            # ADMM-FP
            settings.enable_perturb_admm_fp = True
            settings.enable_restart_admm_fp = True
            settings.k_max_admm_fp_ph1 = 10000 # default
            settings.k_max_admm_fp_ph2 = int(1e7) # large number

            # build and solve motion planning problem
            x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')

            if sol.converged:
                n_feas_admm_fp += 1
                t_admm_fp.append(sol.run_time)
            else:
                t_admm_fp.append(np.nan)

            # ADMM
            settings.enable_perturb_admm_fp = False
            settings.enable_restart_admm_fp = False
            settings.k_max_admm_fp_ph1 = int(1e7) # large number
            settings.k_max_admm_fp_ph2 = 0

            # build and solve motion planning problem
            x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')

            if sol.converged:
                n_feas_admm += 1
                t_admm.append(sol.run_time)
            else:
                t_admm.append(np.nan)

            # run with Gurobi for comparison
            # build and solve motion planning problem
            x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'gurobi')

            if sol.converged:
                n_feas_gurobi += 1
                t_gurobi.append(sol.run_time)
            else:
                t_gurobi.append(np.nan)

            # run with OFP for comparison
            # build and solve motion planning problem
            x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'ofp')

            if sol.converged:
                n_feas_ofp += 1
                t_ofp.append(sol.run_time)
            else:
                t_ofp.append(np.nan)
            

        # log results
        n_feas_admm_fp_arr.append(n_feas_admm_fp)
        n_feas_admm_arr.append(n_feas_admm)
        n_feas_gurobi_arr.append(n_feas_gurobi)
        n_feas_ofp_arr.append(n_feas_ofp)
        t_admm_fp_arr.append(t_admm_fp)
        t_admm_arr.append(t_admm)
        t_gurobi_arr.append(t_gurobi)
        t_ofp_arr.append(t_ofp)

        # print results
        print(f'ADMM-FP feasible in {n_feas_admm_fp} / {n_seeds} cases -> {n_feas_admm_fp/(n_seeds)*100:.2f}%, average solution time = {np.mean(t_admm_fp_arr[-1]) if len(t_admm_fp_arr[-1]) > 0 else np.nan} s')
        print(f'ADMM feasible in {n_feas_admm} / {n_seeds} cases -> {n_feas_admm/(n_seeds)*100:.2f}%, average solution time = {np.mean(t_admm_arr[-1]) if len(t_admm_arr[-1]) > 0 else np.nan} s')
        print(f'Gurobi feasible in {n_feas_gurobi} / {n_seeds} cases -> {n_feas_gurobi/(n_seeds)*100:.2f}%, average solution time = {np.mean(t_gurobi_arr[-1]) if len(t_gurobi_arr[-1]) > 0 else np.nan} s')
        print(f'OFP feasible in {n_feas_ofp} / {n_seeds} cases -> {n_feas_ofp/(n_seeds)*100:.2f}%, average solution time = {np.mean(t_ofp_arr[-1]) if len(t_ofp_arr[-1]) > 0 else np.nan} s')

    # Save results as JSON
    results_dict = {
        'sample_factor_arr': sample_factor_arr,
        'n_seeds': n_seeds,
        'n_feas_admm_fp_arr': n_feas_admm_fp_arr,
        'n_feas_admm_arr': n_feas_admm_arr,
        'n_feas_gurobi_arr': n_feas_gurobi_arr,
        'n_feas_ofp_arr': n_feas_ofp_arr,
        't_admm_fp_arr': t_admm_fp_arr,
        't_admm_arr': t_admm_arr,
        't_gurobi_arr': t_gurobi_arr,
        't_ofp_arr': t_ofp_arr
    }

    with open('reach_avoid_results.json', 'w') as json_file:
        json.dump(results_dict, json_file, indent=4)

    

    ### plots ###
    try:
        subprocess.run(['python', './plot_heuristic_data.py'])
    except Exception as e:
        print(f'Error running plotting script: {e}')

    # generate plots
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

        # figure
        settings.enable_perturb_admm_fp = True
        settings.k_max_admm_fp_ph1 = 10000 # default
        settings.k_max_admm_fp_ph2 = int(1e7) # large number
        sample_factor = 5
        seed_arr = [33, 66, 99]

        figwidth_pt = 505.89
        figsize = (figwidth_pt * inches_per_pt, 0.4*figwidth_pt * inches_per_pt)  # Convert pt to inches
        fig = plt.figure(constrained_layout=True, figsize=figsize)

        width_col_pt = (figwidth_pt / 3.) - textwidth_pt
        gs = fig.add_gridspec(nrows=1, ncols=3, figure=fig, width_ratios=[1,1,1], height_ratios=[1])

        for row in range(1):
            for col in range(3):

                # build and solve motion planning problem
                settings.enable_perturb_admm_fp = True
                seed = seed_arr[row*3 + col]
                print(f'Plotting map for seed {seed}')
                converged = False
                while not converged:
                    x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')
                    converged = sol.converged

                ax = fig.add_subplot(gs[row,col])

                # time-varying state constraints (spatio-temporal corridor)
                zono.plot(Z_map, ax=ax, color='b', alpha=0.25, edgecolor='b', linewidth=1.0)
                zono.plot(zono.intersection(zono.project_onto_dims(XN, [0,1]), Z_map), ax=ax, color='g', alpha=0.5, edgecolor='g', linewidth=1.0)

                # solution
                x_vec = np.array(x_traj)
                ax.plot(x_vec[:,0], x_vec[:,1], '.r', markersize=5)
                
                # axes
                ax.axis('equal')
                ax.grid(alpha=0.2)

                # title
                ax.set_title(f'seed = {seed}', fontsize=textwidth_pt)

        # axis labels
        fig.supxlabel(r'$x$ [m]', fontsize=textwidth_pt)
        fig.supylabel(r'$y$ [m]', fontsize=textwidth_pt)

        # show
        plt.show()

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
    x_traj, u_traj, Z_mpc, Z_map, sol, XN = reach_avoid_problem(seed, sample_factor, settings, 'zonoopt')

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