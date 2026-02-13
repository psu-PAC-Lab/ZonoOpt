import zonoopt as zono
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import subprocess
import copy

# flag indicating what this script does
# options: 'one-shot', 'closed-loop', 'one-shot-comparison', 'closed-loop-comparison'
MODE = 'one-shot'

if MODE == 'one-shot-comparison' or MODE == 'closed-loop-comparison':
    import osqp # used version 0.6.7.post3 in paper
    import gurobipy as gp # used version 12.0.1 in paper
    import shapely
    import pandas as pd


class HrepPoly:
    def __init__(self, A=None, b=None):
        self.A = A
        self.b = b

def vrep_2_hrep(V):
    """Convert V-rep polytope to H-rep polytope"""
    
    # make shapely polygon
    P = shapely.geometry.Polygon(V)

    # make sure vertices oriented counter-clockwise
    if not P.exterior.is_ccw:
        P = shapely.geometry.polygon.orient(P, sign=1.0)

    # get vertices
    V = np.array(P.exterior.coords.xy).transpose()
    V = V[0:-1,:] # remove duplicate vertex

    # get H-rep
    nV = V.shape[0] # number of vertices
    A = np.zeros((nV, 2)) # init
    b = np.zeros(nV) # init
    for i in range(nV):
        xi, yi = V[i]
        xip1, yip1 = V[(i+1) % nV]
        dx = xip1 - xi
        dy = yip1 - yi
        A[i] = [dy, -dx]
        b[i] = dy*xi - dx*yi
        Ai_norm = np.linalg.norm(A[i])
        if Ai_norm > 1e-15:
            A[i] /= Ai_norm
            b[i] /= Ai_norm

    return HrepPoly(A, b)

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

class MPCData:
    def __init__(self, x0, xr, A, B, Q, R, N, X_arr, U, settings=zono.OptSettings(), sol=zono.OptSolution(), O_arr=[]):
        """
        MPC data structure constructor

        Args:
            x0: initial state
            xr: reference trajectory (N+1 x n numpy array)
            A = state transition matrix
            B = control matrix
            Q = state cost matrix
            R = control cost matrix
            N = prediction horizon
            X_arr = array of state constraints
            U = control constraint set
            settings = zonoopt.OptSettings object
            sol = zonoopt.OptSolution object
            O_arr = obstacle array object for plotting

        Returns:
            MPCData object
        """
        self.x0 = x0
        self.xr = xr
        self.A = A
        self.B = B
        self.Q = Q
        self.R = R
        self.N = N
        self.X_arr = X_arr
        self.U = U
        self.settings = settings
        self.sol = sol
        self.O_arr = O_arr


class MPCProb:
    def __init__(self, solver_flag='zonoopt', use_hrep=False, sample_factor=1, verbose=True, prediction_horizon=None):
        """
        Build MPC problem for example

        Args:
            solver_flag (str): options: 'zonoopt', 'osqp', 'gurobi'
            use_hrep (bool): whether to build problem in H-rep form (True) or CZ form (False)
            sample_factor (float): factor to increase/decrease number of samples along path, which affects the fineness of the time-varying constraint sets
            verbose (bool): print information to terminal
            prediction_horizon (int): prediction horizon length (overrides default value)

        Returns:
            data (MPCData): object containing problem parameters and solver settings
        """

        # flags
        self.solver_flag = solver_flag
        self.use_hrep = use_hrep
        self.sample_factor = sample_factor
        self.verbose = verbose
        self.prediction_horizon = prediction_horizon

        ### Environment setup ###

        # path plan
        path_plan = np.array([[0.0, -10.0],
                            [0.0, 15.0],
                            [15.0, 20.0],
                            [15.0, 35.0]])

        # get points at uniform distance traveled along path
        d_path = np.hstack([0.0, np.cumsum(np.sqrt(np.sum(np.diff(path_plan, axis=0)**2, axis=1)))])
        ds = 1.0/self.sample_factor
        d_vec = np.linspace(0.0, d_path[-1], int(d_path[-1]/ds))
        path_x = np.interp(d_vec, d_path, path_plan[:, 0])
        path_y = np.interp(d_vec, d_path, path_plan[:, 1])

        # get indices for the part of the path in the middle
        idx_middle = np.where((path_y >= 15.0) & (path_y <= 20.0))[0]

        # obstacles
        m = (20.0-15.0)/(15.0-0.0) # slope
        V_O1 = np.array([[5.0, 15.6],
                    [10.0, 15.6 + 5.0*m],
                    [4.0, 10.6],
                    [15.0, 12.1],
                    [13.0, 6.6],
                    [7.5, 8.1],
                    [7.2, 16.8]])
        O1 = zono.vrep_2_conzono(V_O1)

        V_02 = np.array([[4.0, 17.6],
                        [9.0, 17.6 + 5.0*m],
                        [1.0, 21.4],
                        [2.0, 26.4],
                        [6.0, 30.4],
                        [11.0, 26.4],
                        [6.9, 17.9]])
        O2 = zono.vrep_2_conzono(V_02)

        # base zonotope for constraint sets
        Z_base = zono.make_regular_zono_2D(1., 6)

        # variable constraint sets
        len_front = int(len(idx_middle)/2)
        len_back = len(idx_middle) - len_front
        rt_vec_middle = np.hstack([np.linspace(2.0, 0.4, len_front),
                            np.linspace(0.4, 2.0, len_back)])
        rs_vec_middle = np.hstack([np.linspace(2.0, 2.0, len_front),
                            np.linspace(2.0, 2.0, len_back)])

        rt_vec = 4.5*np.ones(len(d_vec))
        rs_vec = 4.5*np.ones(len(d_vec))
        rt_vec[idx_middle] = rt_vec_middle
        rs_vec[idx_middle] = rs_vec_middle

        phi_vec = np.hstack(np.arctan2(np.diff(path_y), np.diff(path_x)))
        phi_vec = np.hstack([phi_vec[0], phi_vec])

        Z_cons_arr = []
        for i in range(len(d_vec)):
            C_rot = sparse.csc_matrix([[np.cos(phi_vec[i]), -np.sin(phi_vec[i])],
                                    [np.sin(phi_vec[i]), np.cos(phi_vec[i])]])
            C_lin = sparse.diags([rs_vec[i], rt_vec[i]])
            C = C_rot*C_lin
            b = np.array([path_x[i], path_y[i]])

            Z_cons_arr.append(zono.affine_map(Z_base, C, b))


        ### MPC setup ###

        # time step
        dt = 1.0/self.sample_factor

        # 2D double integrator dynamics
        # x = [x, y, xdot, y_dot]
        # u = [x_ddot, y_ddot]
        A = sparse.csc_matrix(
            [[1., 0., dt, 0.],
            [0., 1., 0., dt],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.]])
        B = sparse.csc_matrix(
            [[0.5*dt**2, 0.],
            [0., 0.5*dt**2],
            [dt, 0.],
            [0., dt]])

        # state feasible set

        # velocity constraints
        v_max = 5.0
        Xv = zono.make_regular_zono_2D(v_max, 12)

        # input feasible set

        # turn rate constraints
        om_max = 75*np.pi/180
        v_min = 0.1 # fictitious, om_max not enforced below this
        U = zono.make_regular_zono_2D(om_max*v_min, 12)

        # MPC horizon
        N = 55*self.sample_factor

        # cost function matrices
        Q = sparse.diags([1.0, 1.0, 0.0, 0.0])
        R = 10*sparse.eye(2)

        # solver settings
        settings = zono.OptSettings()
        settings.inf_norm_conv = True
        settings.rho = 1.
        settings.eps_prim = 1e-2
        settings.eps_dual = 1e-2

        ### Solve example ###

        # sim length
        n_sim = len(d_vec)

        # reference
        x_ref = np.zeros((0, 4)) # init
        for i in range(n_sim):
            x_ref = np.vstack((x_ref, np.array([path_x[i], path_y[i], 0.0, 0.0])))
        for k in range(N+1):
            x_ref = np.vstack((x_ref, x_ref[-1,:]))
        x_ref = x_ref.transpose()

        # initial condition
        x = np.array([0.0, -10.0, 0.0, 0.0]) # init

        # reference
        xr = []
        for i in range(1,N+1):
            xr.append(x_ref[:,i].flatten())

        # build constraint sets
        X_arr = []
        for i in range(N+1):
            if i >= len(Z_cons_arr):
                X_arr.append(zono.cartesian_product(Z_cons_arr[-1], Xv))
            else:
                X_arr.append(zono.cartesian_product(Z_cons_arr[i], Xv))

        # handle H-rep case
        if self.use_hrep:

            # build constraint sets in H-rep
            X_arr_hrep = []
            for i in range(N+1):
                Xv_hrep = vrep_2_hrep(zono.get_vertices(Xv))
                if i >= len(Z_cons_arr):
                    Xp_hrep = vrep_2_hrep(zono.get_vertices(Z_cons_arr[-1]))
                else:
                    Xp_hrep = vrep_2_hrep(zono.get_vertices(Z_cons_arr[i]))
                A_hrep = sparse.block_diag((sparse.csc_matrix(Xp_hrep.A), sparse.csc_matrix(Xv_hrep.A)))
                b_hrep = np.hstack((Xp_hrep.b, Xv_hrep.b))
                X_arr_hrep.append(HrepPoly(A_hrep, b_hrep))

            U_hrep = vrep_2_hrep(zono.get_vertices(U))

            self.data = MPCData(x, xr, A, B, Q, R, N, X_arr_hrep, U_hrep, settings=settings, O_arr=(O1, O2))

        else:
            self.data = MPCData(x, xr, A, B, Q, R, N, X_arr, U, settings=settings, O_arr=(O1, O2))

        # initial condition
        self.x = x

    def simulate(self):
        """
        Simulate closed-loop system with MPC controller

        Returns:
            x_sim (list of numpy arrays): state trajectory
            u_sim (list of numpy arrays): control trajectory
            x_traj_arr (list of lists of numpy arrays): MPC state trajectories at each time step
            u_traj_arr (list of lists of numpy arrays): MPC control trajectories at each time step
            sol_time_arr (list of floats): MPC solution times at each time step
            iter_arr (list of ints): number of iterations taken by solver at each time step
        """

        # number of simulation time steps
        n_sim = self.data.N

        # initial mpc data structure
        data = copy.copy(self.data)
        if self.prediction_horizon is None:
            data.N = self.data.N
        elif self.prediction_horizon < self.data.N:
            data.N = self.prediction_horizon
        else:
            raise ValueError('Prediction horizon cannot be greater than total problem horizon')
        self._update_mpc_data(data, 0)

        # storage
        x_sim = [self.x]
        u_sim = []
        x_traj_arr = []
        u_traj_arr = []
        self.sim_sol_struct = []

        # simulate
        for k in range(n_sim):
            # solve mpc
            x_traj, u_traj = self.solve_mpc(data)

            # apply first control input
            u_k = u_traj[0]
            self.x = self.data.A.dot(self.x) + self.data.B.dot(u_k)

            # store
            x_sim.append(self.x)
            u_sim.append(u_k)
            x_traj_arr.append(x_traj)
            u_traj_arr.append(u_traj)

            # update mpc data
            self._update_mpc_data(data, k+1)

            # store solution struct
            self.sim_sol_struct.append(data.sol.copy())

        sol_time_arr = np.array([sol.run_time for sol in self.sim_sol_struct])
        iter_arr = np.array([sol.iter for sol in self.sim_sol_struct])

        # log
        if self.verbose:
            print('Closed-loop simulation complete')

            # statistics
            print(f'Average solution time: {np.mean(sol_time_arr)} s, min = {np.min(sol_time_arr)} s, max = {np.max(sol_time_arr)} s')
            print(f'Average iterations: {np.mean(iter_arr)}, min = {np.min(iter_arr)}, max = {np.max(iter_arr)}')

        return x_sim, u_sim, sol_time_arr, iter_arr, x_traj_arr, u_traj_arr

    def solve_mpc(self, data=None):
        """
        Solves MPC problem

        Args:
            data (MPCData): MPC problem data

        Returns:
            x_traj (list of numpy arrays): state trajectory
            u_traj (list of numpy arrays): control trajectory
        """

        if data is None:
            data = self.data

        # solve MPC
        if self.use_hrep:
            x_traj, u_traj = self._hrep_mpc(data)
        else:
            x_traj, u_traj = self._zono_mpc(data)

        # display
        if self.verbose:
            print(f'iter: {data.sol.iter}, run time: {data.sol.run_time}, startup time = {data.sol.startup_time}')

        return x_traj, u_traj

    def plot_traj(self, x_traj):
        """
        Plot trajectory

        Args:
            x_traj (list of numpy arrays): state trajectory
        """

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
        figsize = (245.71 * inches_per_pt, 0.7*245.71 * inches_per_pt)  # Convert pt to inches

        # plot
        with plt.rc_context(rc_context):

            # flip x and y axes
            C_flip = sparse.csc_matrix([[0, 1],
                                [1, 0]])

            # figure
            fig = plt.figure(constrained_layout=True, figsize=figsize)
            ax = fig.add_subplot(111)

            # path plan
            xr_x = [xr_k[0] for xr_k in self.data.xr]
            xr_y = [xr_k[1] for xr_k in self.data.xr]
            ax.plot(xr_y, xr_x, '-k')

            # obstacles
            for O in self.data.O_arr:
                zono.plot(zono.affine_map(O, C_flip), ax=ax, color='g', alpha=0.7, edgecolor='k')

            # time-varying state constraints (spatio-temporal corridor)
            for X in self.data.X_arr:
                zono.plot(zono.affine_map(zono.project_onto_dims(X, [0,1]), C_flip), ax=ax, color='b', alpha=0.05, edgecolor='b', linewidth=1.0)

            # MPC solution
            x_vec = np.array(x_traj)
            ax.plot(x_vec[:,1], x_vec[:,0], '.r', markersize=3)

            # axes
            ax.set_xlabel('$x$ [m]')
            ax.set_ylabel('$y$ [m]')
            ax.axis('equal')
            ax.grid(alpha=0.2)

            # save
            if is_latex_installed():
                plt.savefig('mpc_time_varying_cons.pgf')

            plt.show()

    def _update_mpc_data(self, data, k):
        """
        Update MPC data structure for time step k

        Args:
            data (MPCData): current MPC data structure
            k (int): current time step
        """

        # update initial condition
        data.x0 = self.x

        # update reference trajectory
        data.xr = []
        for i in range(1, data.N+1):
            if (k+i) < len(self.data.xr):
                data.xr.append(self.data.xr[k+i])
            else:
                data.xr.append(self.data.xr[-1])

        # update constraint sets
        data.X_arr = []
        for i in range(data.N+1):
            if (k+i) < len(self.data.X_arr):
                data.X_arr.append(self.data.X_arr[k+i])
            else:
                data.X_arr.append(self.data.X_arr[-1])

    def _zono_mpc(self, data):
        """
        Build and solve MPC problem in CZ form.

        Args:
            data (MPCData): MPC problem data

        Returns:
            x_traj (list of numpy arrays): state trajectory
            u_traj (list of numpy arrays): control trajectory
            Z (zonoopt.ConZono): final constraint set after MPC formulation
        """

        # dims
        nx = data.A.shape[1]
        nu = data.B.shape[1]

        # index tracking
        idx = 0
        idx_x = []
        idx_u = []

        # initial state
        Z = zono.Point(data.x0)
        idx_x.append([j for j in range(idx, idx+nx)])
        idx += nx

        # init cost
        P = data.Q
        q = np.zeros(nx)

        for k in range(data.N):

            # control
            Z = zono.cartesian_product(Z, data.U)

            idx_u.append([j for j in range(idx, idx+nu)])
            idx += nu

            # dynamics
            nZ = Z.get_n()
            Z = zono.cartesian_product(Z, data.X_arr[k+1])
            ABmI = sparse.hstack((sparse.csc_matrix((nx, nZ-(nx+nu))), data.A, data.B, -sparse.eye(nx)))
            Z = zono.intersection(Z, zono.Point(np.zeros(nx)), ABmI)

            idx_x.append([j for j in range(idx, idx+nx)])
            idx += nx

            # cost
            P = sparse.block_diag((P, data.R, data.Q))
            q = np.hstack([q, np.zeros(nu), -data.Q.dot(data.xr[k])])

        if self.solver_flag=='zonoopt':

            # optimize
            xopt = Z.optimize_over(P, q, settings=data.settings, solution=data.sol)

        elif self.solver_flag=='osqp':

            # build osqp problem
            P_tilde = Z.get_G().transpose().dot(P.dot(Z.get_G()))
            q_tilde = Z.get_G().transpose().dot(P*Z.get_c() + q)
            A_osqp = Z.get_A()
            A_osqp = sparse.vstack([A_osqp, sparse.eye(Z.get_nG())])
            if Z.is_0_1_form():
                l_osqp = np.hstack([Z.get_b(), np.zeros(Z.get_nG())])
            else:
                l_osqp = np.hstack([Z.get_b(), -1*np.ones(Z.get_nG())])
            u_osqp = np.hstack([Z.get_b(), np.ones(Z.get_nG())])

            # create object
            prob = osqp.OSQP()
            prob.setup(P=P_tilde, q=q_tilde, A=A_osqp, l=l_osqp, u=u_osqp, eps_rel=0.0, eps_abs=data.settings.eps_prim)

            # solve
            osqp_sol = prob.solve()
            xopt = Z.get_G()*osqp_sol.x + Z.get_c()

            # logging
            data.sol.run_time = osqp_sol.info.run_time
            data.sol.iter = osqp_sol.info.iter
            data.sol.startup_time = osqp_sol.info.setup_time
            data.sol.infeasible = osqp_sol.info.status_val != 1

        elif self.solver_flag=='gurobi':

            # build gurobi model
            P_tilde = Z.get_G().transpose().dot(P.dot(Z.get_G()))
            q_tilde = Z.get_G().transpose().dot(P*Z.get_c() + q)

            # create model
            prob = gp.Model()

            # add variables
            if Z.is_0_1_form():
                x_gurobi = prob.addMVar(Z.get_nG(), lb=np.zeros(Z.get_nG()), ub=np.ones(Z.get_nG()))
            else:
                x_gurobi = prob.addMVar(Z.get_nG(), lb=-np.ones(Z.get_nG()), ub=np.ones(Z.get_nG()))

            # add constraints
            A_gurobi = Z.get_A()
            b_gurobi = Z.get_b()
            prob.addConstr(A_gurobi.dot(x_gurobi) == b_gurobi)

            # add objective
            prob.setMObjective(P_tilde, q_tilde, 0.0, sense=gp.GRB.MINIMIZE)

            # optimize
            prob.optimize()
            xi_opt = np.array([x.X for x in prob.getVars()])
            xopt = Z.get_G()*xi_opt + Z.get_c()

            # logging
            data.sol.run_time = prob.Runtime
            data.sol.iter = prob.BarIterCount
            data.sol.startup_time = 0.0
            data.sol.infeasible = prob.Status != 2

        else:
            return ValueError('Unknown solver flag')

        if self.verbose:
            print(f'Z.n: {Z.get_n()}, Z.nG: {Z.get_nG()}, Z.nC: {Z.get_nC()}')
            print(f'NNZ(G): {Z.get_G().nnz}, NNZ(A): {Z.get_A().nnz}')

        # state, input trajectory
        x_traj = []
        u_traj = []
        for idx in idx_x:
            x_traj.append(xopt[idx])
        for idx in idx_u:
            u_traj.append(xopt[idx])

        # return trajectories
        return x_traj, u_traj

    def _hrep_mpc(self, data):
        """
        Build and solve MPC problem in H-rep form.

        Args:
            data (MPCData): MPC problem data

        Returns:
            x_traj (list of numpy arrays): state trajectory
            u_traj (list of numpy arrays): control trajectory
        """

        # problem dimensions
        nx = data.A.shape[1]
        nu = data.B.shape[1]

        # index tracking
        idx = 0
        idx_x = []
        idx_u = []

        # build equality constraint matrix
        C = -sparse.eye(nx)
        idx_x.append([j for j in range(idx, idx+nx)])
        idx += nx

        for k in range(data.N):
            C = sparse.hstack((C, sparse.csc_matrix((C.shape[0], nu))))
            AB = sparse.hstack((sparse.csc_matrix((nx, C.shape[1]-nx-nu)), sparse.hstack((data.A, data.B))))
            C = sparse.vstack((C, AB))

            idx_u.append([j for j in range(idx, idx+nu)])
            idx += nu

            mI = sparse.vstack((sparse.csc_matrix((C.shape[0]-nx, nx)), -sparse.eye(nx)))
            C = sparse.hstack((C, mI))

            idx_x.append([j for j in range(idx, idx+nx)])
            idx += nx

        # equality constraint vector
        d = np.zeros(C.shape[0])
        d[0:nx] = -data.x0

        # build inequality constraint matrix
        G = sparse.hstack((sparse.csc_matrix((data.U.A.shape[0],nx)), data.U.A))
        for k in range(1,data.N+1):
            if k < data.N:
                G = sparse.block_diag((G, data.X_arr[k].A, data.U.A))
            else:
                G = sparse.block_diag((G, data.X_arr[k].A))

        # inequality constraint vector
        w = data.U.b
        for k in range(1, data.N+1):
            if k < data.N:
                w = np.hstack((w, data.X_arr[k].b, data.U.b))
            else:
                w = np.hstack((w, data.X_arr[k].b))

        # cost function
        P = sparse.csc_matrix((0,0))
        for k in range(data.N+1):
            if k < data.N:
                P = sparse.block_diag((P, data.Q, data.R))
            else:
                P = sparse.block_diag((P, data.Q))

        q = np.zeros(nx+nu)
        for k in range(1, data.N+1):
            if k < data.N:
                q = np.hstack((q, -data.Q.dot(data.xr[k-1]), np.zeros(nu)))
            else:
                q = np.hstack((q, -data.Q.dot(data.xr[k-1])))

        # ensure all matrices are in CSC format
        C = C.tocsc()
        G = G.tocsc()
        P = P.tocsc()

        # solve
        if self.solver_flag=='zonoopt':
            raise ValueError('Not implemented')

        elif self.solver_flag=='osqp':

            # build osqp problem
            A_osqp = C
            A_osqp = sparse.vstack([A_osqp, G])
            l_osqp = np.hstack([d, -np.inf*np.ones(len(w))])
            u_osqp = np.hstack([d, w])

            # create object
            prob = osqp.OSQP()
            prob.setup(P=P, q=q, A=A_osqp, l=l_osqp, u=u_osqp, eps_rel=0.0, eps_abs=data.settings.eps_prim)

            # solve
            osqp_sol = prob.solve()
            xopt = osqp_sol.x

            # logging
            data.sol.run_time = osqp_sol.info.run_time
            data.sol.iter = osqp_sol.info.iter
            data.sol.startup_time = osqp_sol.info.setup_time
            data.sol.infeasible = osqp_sol.info.status_val != 1

        elif self.solver_flag=='gurobi':

            # create model
            prob = gp.Model()

            # add variables
            n_vars = len(q)
            x_gurobi = prob.addMVar(n_vars, lb=-gp.GRB.INFINITY*np.ones(n_vars), ub=gp.GRB.INFINITY*np.ones(n_vars))

            # add constraints
            prob.addConstr(C.dot(x_gurobi) == d)
            prob.addConstr(G.dot(x_gurobi) <= w)

            # add objective
            prob.setMObjective(P, q, 0.0, sense=gp.GRB.MINIMIZE)

            # optimize
            prob.optimize()
            xopt = np.array([x.X for x in prob.getVars()])

            # logging
            data.sol.run_time = prob.Runtime
            data.sol.iter = prob.BarIterCount
            data.sol.startup_time = 0.0
            data.sol.infeasible = prob.Status != 2

        else:
            raise ValueError('Unknown solver flag')

        # return state, input trajectory
        x_traj = []
        u_traj = []
        for idx in idx_x:
            x_traj.append(xopt[idx])
        for idx in idx_u:
            u_traj.append(xopt[idx])

        if self.verbose:
            print(f'P dims = {P.shape}, C dims = {C.shape}, G dims = {G.shape}')
            print(f'NNZ(P): {P.nnz}, NNZ(C): {C.nnz}, NNZ(G): {G.nnz}')

        # return trajectories
        return x_traj, u_traj

    

### run example ###
if  MODE == 'one-shot':

    # build problem
    prob = MPCProb()

    # solve
    x_traj_sol, u_traj_sol = prob.solve_mpc()

    # plot solution
    prob.plot_traj(x_traj_sol)

elif MODE == 'closed-loop':

    # build problem
    prob = MPCProb(sample_factor=4, prediction_horizon=100, solver_flag='zonoopt')

    # solve
    x_traj_sol, u_traj_sol, sol_time_arr, _, x_traj_arr, _ = prob.simulate()

    # plot simulated trajectory
    prob.plot_traj(x_traj_sol) # passing x_traj_arr[k] for a given k plots the MPC solution at that time step

    # plot solution times
    plt.figure()
    ax = plt.subplot(111)
    ax.plot(sol_time_arr*1000., '-o')
    ax.set_xlabel('Time step k')
    ax.set_ylabel('[ms]')
    ax.set_title('MPC solution time vs. time step')
    ax.grid(alpha=0.2)
    plt.show()

elif MODE == 'one-shot-comparison':

    # sample factors
    sample_factor_arr = []
    for i in range(1,21+1):
        sample_factor_arr.append(i)

    # solver options: 'zonoopt', 'osqp', 'gurobi'
    solver_flag_arr = ['zonoopt', 'osqp', 'gurobi']

    # representation options
    use_hrep_arr = [False, True]

    # solution struct logging
    sol_log = pd.DataFrame()

    # run example in loop
    for sample_factor in sample_factor_arr:
        for solver_flag in solver_flag_arr:
            for use_hrep in use_hrep_arr:

                # solve
                try:
                    prob = MPCProb(solver_flag=solver_flag, use_hrep=use_hrep, sample_factor=sample_factor, verbose=False)
                    prob.solve_mpc()
                except ValueError as e:
                    continue

                # log solution
                new_element = pd.DataFrame({'sample_factor': sample_factor,
                                            'self.solver_flag': solver_flag,
                                            'self.use_hrep': use_hrep,
                                            'sol': prob.data.sol.copy()}, index=[0])
                sol_log = pd.concat([sol_log, new_element], ignore_index=True)

    # plot solution times
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
    figsize = (245.71 * inches_per_pt, 0.9*245.71 * inches_per_pt)  # Convert pt to inches

    # plot
    with plt.rc_context(rc_context):

        fig = plt.figure(constrained_layout=True, figsize=figsize)
        ax = fig.add_subplot(111)

        # get solution time vs sample factor for each solver
        for solver_flag in solver_flag_arr:
            for use_hrep in use_hrep_arr:

                # non-existent case
                if solver_flag == 'zonoopt' and use_hrep:
                    continue

                # label
                if solver_flag == 'zonoopt':
                    solver_str = r'ZonoOpt'
                elif solver_flag == 'osqp':
                    solver_str = r'OSQP'
                elif solver_flag == 'gurobi':
                    solver_str = r'Gurobi'
                
                if use_hrep:
                    method_str = r'H-rep'
                else:
                    method_str = r'CZ'

                label = solver_str + r' ' + method_str

                # filter
                sol_log_filt = sol_log[(sol_log['solver_flag'] == solver_flag) & (sol_log['use_hrep'] == use_hrep)]

                # plot
                ax.plot(sol_log_filt['sample_factor']*55, 1000.*sol_log_filt['sol'].apply(lambda x: x.run_time), label=label)

        # labels
        ax.set_xlabel(r'Horizon $N$')
        ax.set_ylabel(r'Solution time [ms]')
        ax.legend()

        # x ticks
        xticks = [55*i for i in range(1,21+1)]
        xticklabels = []
        for i in range(1,21+1):
            if (i-1) % 5 == 0:
                xticklabels.append(str(55*i))
            else:
                xticklabels.append('')
        ax.set_xticks(xticks, xticklabels)

        # grid on
        ax.grid(alpha=0.2)

        # save
        if is_latex_installed():
            plt.savefig('solution_method_comparison.pgf')

        # display
        plt.show()

elif MODE == 'closed-loop-comparison':

    cases = ('ZonoOpt CZ', 'OSQP CZ', 'OSQP H-rep', 'Gurobi CZ', 'Gurobi H-rep')

    sol_times_ms = []
    for case in cases:

        # build problem
        if case == 'ZonoOpt CZ':
            solver_flag = 'zonoopt'
            use_hrep = False
        elif case == 'OSQP CZ':
            solver_flag = 'osqp'
            use_hrep = False
        elif case == 'OSQP H-rep':
            solver_flag = 'osqp'
            use_hrep = True
        elif case == 'Gurobi CZ':
            solver_flag = 'gurobi'
            use_hrep = False
        elif case == 'Gurobi H-rep':
            solver_flag = 'gurobi'
            use_hrep = True
        else:
            raise ValueError('Unknown case')
        prob = MPCProb(sample_factor=4, prediction_horizon=100, solver_flag=solver_flag, use_hrep=use_hrep, verbose=False)

        # solve
        x_traj_sol, u_traj_sol, sol_time_arr, iter_arr, _, _ = prob.simulate()
        sol_times_ms.append(1000. * sol_time_arr) # convert to ms

    # box plot of solution times
    textwidth_pt = 10.
    rc_context = {
        "text.usetex": True,
        "font.size": textwidth_pt,
        "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
        "pgf.texsystem": "pdflatex",
        "pgf.rcfonts": False,
    }
    inches_per_pt = 1 / 72.27
    figsize = (245.71 * inches_per_pt, 0.7*245.71 * inches_per_pt)  # Convert pt to inches

    with plt.rc_context(rc_context):

        fig = plt.figure(constrained_layout=True, figsize=figsize)
        ax = fig.add_subplot(111)

        labels = (r'\begin{tabular}{c} ZonoOpt \\ CZ \end{tabular}', r'\begin{tabular}{c} OSQP \\ CZ \end{tabular}', r'\begin{tabular}{c} OSQP \\ H-rep \end{tabular}', r'\begin{tabular}{c} Gurobi \\ CZ \end{tabular}', r'\begin{tabular}{c} Gurobi \\ H-rep \end{tabular}')

        bp = ax.boxplot(
            sol_times_ms,
            tick_labels=labels,
            patch_artist=True,
            medianprops={'color': 'black'},
            showfliers=True,
            sym='k.',
            whis=(0., 100.) # whiskers cover all data
        )
        ax.set_ylabel(r'Solution time [ms]')
        ax.set_yscale('log')
        ax.set_ylim((1, 100))
        ax.grid(alpha=0.2)

        colors = ['blue', 'orange', 'green', 'red', 'purple']
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i % len(colors)])
            box.set_alpha(0.5)


        # save
        if is_latex_installed():
            plt.savefig('closed_loop_solution_time_comparison.pgf')

        plt.show()



else:
    raise ValueError('Unknown mode flag')