import zonoopt as zono
from scipy import sparse
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

# Ball and paddle system from Marcucci, Tobia, and Russ Tedrake. "Mixed-integer formulations for optimal control of piecewise-affine systems." Proceedings of the 22nd ACM International Conference on Hybrid Systems: Computation and Control. 2019.
# Modes and associated domains generated in https://github.com/jrobbins11/pympc/blob/feature/bouncing-ball-zono/repeatability_evaluation/pwa_dynamics.py
# - this is a fork of https://github.com/TobiaMarcucci/pympc


# constants
l = 0.6 # paddle / ceiling length
d = 0.4 # ceiling height
r = 0.1 # ball radius
dt = 0.1 # time step

# mode
MODE = 'representative' # 'multi-run' or 'representative'

class AffineSystem:
    def __init__(self, A, B, c):
        self.A = A
        self.B = B
        self.c = c

    def to_json(self, filename):
        dict = {
            "A": AffineSystem._sparse_to_dict(sparse.coo_matrix(self.A)),
            "B": AffineSystem._sparse_to_dict(sparse.coo_matrix(self.B)),
            "c": self.c.tolist()
        }
        with open(filename, 'w') as f:
            json.dump(dict, f)

    @staticmethod
    def from_json(filename):
        with open(filename, 'r') as f:
            dict = json.load(f)
        A = AffineSystem._dict_to_sparse(dict["A"]).toarray()
        B = AffineSystem._dict_to_sparse(dict["B"]).toarray()
        c = np.array(dict["c"])
        return AffineSystem(A, B, c)

    @staticmethod
    def _sparse_to_dict(m):
        coo = m.tocoo()
        return {
            "row": coo.row.tolist(),
            "col": coo.col.tolist(),
            "data": coo.data.tolist(),
            "shape": coo.shape
        }

    @staticmethod
    def _dict_to_sparse(d):
        return sparse.coo_matrix((d["data"], (d["row"], d["col"])), shape=d["shape"]).tocsc()
    

def make_graph_of_function(A, B, E, SUin, Sout=None):
    """
    Make graph of function Psi for function x+ = A x + B u + E
    
    Args:
        A (np.ndarray): state matrix
        B (np.ndarray): input matrix
        E (np.ndarray): offset vector
        SUin (zono.HybZono): input set for (x, u)
        Sout (zono.HybZono): output set for x+

    Returns:
        zono.HybZono: graph of function Psi
    """

    nx = A.shape[0]
    nu = B.shape[1]
      
    if Sout is None:
        IIAB = np.vstack( (np.eye(nx+nu), np.hstack( (A, B) )))
        E_ext = np.hstack( (np.zeros(nx+nu), E) )
        Psi = zono.minkowski_sum(zono.affine_map(SUin, IIAB), zono.Point(E_ext))
    else:
        Psi = zono.cartesian_product(SUin, Sout)
        Psi = zono.intersection(Psi, zono.Point(-E), np.hstack((A, B, -np.eye(nx))))

    return Psi

    
def zono_planning_prob(x0, xr, Psi_arr, S, U, Q, R, QN, N, settings=zono.OptSettings(), sol=zono.OptSolution()):
    """
    Builds planning problem using graph of function Psi.

    Args:
        x0 (np.ndarray): initial state
        xr (list of np.ndarray): reference state trajectory
        Psi_arr (tuple[zono.HybZono]): graph of function
        S (zono.HybZono): state constraint set
        U (zono.HybZono): input constraint set
        Q (scipy.sparse.csc_matrix): state cost matrix
        R (scipy.sparse.csc_matrix): input cost matrix
        QN (scipy.sparse.csc_matrix): terminal state cost matrix
        N (int): horizon length
        settings (zono.OptSettings, optional): optimization settings
        sol (zono.OptSolution, optional): optimization solution reference

    Returns:
        tuple: (x_traj, u_traj, Z)
            x_traj (list of np.ndarray): state trajectory
            u_traj (list of np.ndarray): input trajectory
            Z (zono.HybZono): optimization problem set
            z (np.ndarray): optimal solution
    """

    # dims
    nx = S.get_n()
    nu = U.get_n()

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

        # lift
        Z = zono.cartesian_product(Z, U)
        Z = zono.cartesian_product(Z, S)

        idx_u.append([j for j in range(idx, idx+nu)])
        idx += nu

        idx_x.append([j for j in range(idx, idx+nx)])
        idx += nx
        
        # dynamics and constraints
        nZ = Z.get_n()
        Z = zono.intersection_over_dims(Z, Psi_arr[k], [i for i in range(nZ - Psi_arr[k].get_n(), nZ)])

        # cost
        if k == N-1:
            P = sparse.block_diag((P, R, QN))
            q = np.hstack([q, np.zeros(nu), -QN.dot(xr[k])])
            c += 0.5*xr[k].dot(QN.dot(xr[k]))
        else:
            P = sparse.block_diag((P, R, Q))
            q = np.hstack([q, np.zeros(nu), -Q.dot(xr[k])])
            c += 0.5*xr[k].dot(Q.dot(xr[k]))

    # solve optimization problem
    xopt = Z.optimize_over(P, q, c=c, settings=settings, solution=sol)

    # state, input trajectory
    x_traj = []
    u_traj = []
    for idx in idx_x:
        x_traj.append(xopt[idx])
    for idx in idx_u:
        u_traj.append(xopt[idx])

    # return
    return x_traj, u_traj, xopt, Z, P, q, c, idx_x, idx_u

def simulate_pwa(x0, u_traj, sys_arr, D_arr):
    """
    Simulate PWA system given mode dynamics and domains

    Args:
        x0 (np.ndarray): initial state
        u_traj (list of np.ndarray): input trajectory
        sys_arr (list of AffineSystem): mode dynamics
        D_arr (list of zono.HybZono): mode domains

    Returns:
        list of np.ndarray: state trajectory
    """

    x_traj = [x0]
    x_k = x0.copy()

    for u_k in u_traj:
        # find mode
        mode = None
        for i in range(len(D_arr)):
            D_i = D_arr[i]
            if D_i.contains_point(np.hstack((x_k, u_k))):
                if mode is None:
                    mode = i
                else:
                    raise ValueError("Multiple modes valid for state and input.")
        if mode is not None:
            sys_i = sys_arr[mode]
            x_kp1 = sys_i.A.dot(x_k) + sys_i.B.dot(u_k) + sys_i.c
            x_traj.append(x_kp1)
            x_k = x_kp1   
        else:
            raise ValueError("No valid mode found for state and input.")

    return x_traj


class BallPlot:

    def __init__(self, x_traj, x_traj_shadow=None):
        self.x_traj = x_traj
        self.x_traj_shadow = x_traj_shadow
        self.x_min, self.x_max, self.y_min, self.y_max = self.get_fig_bounds()
        self.shadow_alpha = 0.3

    def make_ball_paddle_fig(self, ncols=4, k_list=None):
        """
        Make ball and paddle figure
        """

        if k_list is None:
            k_list = range(len(self.x_traj))

        # latex settings
        textwidth_pt = 10.
        rc_context = {
            "text.usetex": True,
            "font.size": textwidth_pt,
            "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
            "pgf.texsystem": "pdflatex",
            "pgf.rcfonts": False,
        }

        inches_per_pt = 1 / 72.27

        with plt.rc_context(rc_context):

            # make figure
            figwidth_pt = 245.71
            figsize = (figwidth_pt * inches_per_pt, 0.8*figwidth_pt * inches_per_pt)  # Convert pt to inches
            fig = plt.figure(constrained_layout=True, figsize=figsize)
            
            nrows = int(np.ceil(len(k_list) / ncols))
            gs = fig.add_gridspec(ncols=ncols, nrows=nrows)

            i = 0
            j = 0
            for k in k_list:
                if j >= ncols:
                    j = 0
                    i += 1

                x_k = x_traj[k]

                ax = fig.add_subplot(gs[i,j])

                BallPlot.draw_ceiling(ax)
                BallPlot.draw_paddle(ax, x_k[3], x_k[4])
                BallPlot.draw_ball(ax, x_k[0], x_k[1], x_k[2])
                if self.x_traj_shadow is not None:
                    x_k_shadow = self.x_traj_shadow[k]
                    BallPlot.draw_ball(ax, x_k_shadow[0], x_k_shadow[1], x_k_shadow[2], alpha=self.shadow_alpha)

                ax.set_aspect('equal')
                ax.set_xlim((self.x_min, self.x_max))
                ax.set_ylim((self.y_min, self.y_max))
                ax.grid(alpha=0.2)

                if j != 0:
                    ax.set_yticklabels([])
                ax.set_xticks([-0.25, 0., 0.25])
                ax.set_xticklabels(['-0.25', '0', '0.25'])
                ax.set_title(r'$k=' + str(k) + '$', fontsize=textwidth_pt)

                j += 1

            fig.supxlabel(r'$x$ [m]', fontsize=textwidth_pt)
            fig.supylabel(r'$y$ [m]', fontsize=textwidth_pt)

            plt.savefig('ball_paddle_traj.pgf', dpi=300)
            plt.show()
            

    def make_ball_paddle_animation(self, filename='ball_paddle.mp4', playback_speed=1.0, fps=30, dpi=300):
        """
        Make animation for ball and paddle example
        """

        # get trajectories for interpolation
        x_ball = [x_k[0] for x_k in self.x_traj]
        y_ball = [x_k[1] for x_k in self.x_traj]
        th_ball = [x_k[2] for x_k in self.x_traj]
        x_paddle = [x_k[3] for x_k in self.x_traj]
        y_paddle = [x_k[4] for x_k in self.x_traj]
        t_arr = np.arange(len(self.x_traj)) * dt

        if self.x_traj_shadow is not None:
            x_ball_shadow = [x_k[0] for x_k in self.x_traj_shadow]
            y_ball_shadow = [x_k[1] for x_k in self.x_traj_shadow]
            th_ball_shadow = [x_k[2] for x_k in self.x_traj_shadow]
            x_paddle_shadow = [x_k[3] for x_k in self.x_traj_shadow]
            y_paddle_shadow = [x_k[4] for x_k in self.x_traj_shadow]
        else:
            x_ball_shadow = None
            y_ball_shadow = None
            th_ball_shadow = None
            x_paddle_shadow = None
            y_paddle_shadow = None

        # latex settings
        textwidth_pt = 10.
        rc_context = {
            "text.usetex": True,
            "font.size": textwidth_pt,
            "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
            "pgf.texsystem": "pdflatex",
            "pgf.rcfonts": False,
        }

        inches_per_pt = 1 / 72.27
        figwidth_pt = 245.71
        figheight_pt = figwidth_pt

        with plt.rc_context(rc_context):

            # make figure
            fig = plt.figure(figsize=(figwidth_pt * inches_per_pt, figheight_pt * inches_per_pt), constrained_layout=True)
            ax = fig.add_subplot(1,1,1)
            ax.set_aspect('equal')

            # animation
            writer = FFMpegWriter(fps=fps)
            writer.setup(fig, filename, dpi=dpi)

            # frames
            t_max_anim = (len(self.x_traj)-1) * dt * playback_speed
            t = 0.
            while t < t_max_anim:
                ax.cla()

                t += 1./fps
                t_interp = t / playback_speed

                x_ball_t = np.interp(t_interp, t_arr, x_ball)
                y_ball_t = np.interp(t_interp, t_arr, y_ball)
                th_ball_t = np.interp(t_interp, t_arr, th_ball)
                x_paddle_t = np.interp(t_interp, t_arr, x_paddle)
                y_paddle_t = np.interp(t_interp, t_arr, y_paddle)

                # axes
                ax.set_xlabel(r'$x$ [m]')
                ax.set_ylabel(r'$y$ [m]')
                ax.set_xlim((self.x_min, self.x_max))
                ax.set_ylim((self.y_min, self.y_max))
                ax.grid(alpha=0.2)

                # ball
                BallPlot.draw_ceiling(ax)
                BallPlot.draw_paddle(ax, x_paddle_t, y_paddle_t)
                BallPlot.draw_ball(ax, x_ball_t, y_ball_t, th_ball_t)

                # shadow
                if self.x_traj_shadow is not None:
                    x_ball_t_shadow = np.interp(t_interp, t_arr, x_ball_shadow)
                    y_ball_t_shadow = np.interp(t_interp, t_arr, y_ball_shadow)
                    th_ball_t_shadow = np.interp(t_interp, t_arr, th_ball_shadow)
                    x_paddle_t_shadow = np.interp(t_interp, t_arr, x_paddle_shadow)
                    y_paddle_t_shadow = np.interp(t_interp, t_arr, y_paddle_shadow)
                    BallPlot.draw_paddle(ax, x_paddle_t_shadow, y_paddle_t_shadow)
                    BallPlot.draw_ball(ax, x_ball_t_shadow, y_ball_t_shadow, th_ball_t_shadow, alpha=self.shadow_alpha)

                fig.canvas.flush_events()
                writer.grab_frame()

    @staticmethod
    def draw_ball(ax, x, y, th, curve_pts=100, alpha=1.0):
        """
        Draw ball for bouncing ball example
        """

        th = -th

        # green section
        th_grn = np.linspace(th, th+np.pi, curve_pts)
        x_grn = x + r*np.sin(th_grn)
        y_grn = y + r + r*np.cos(th_grn)

        # blue section
        th_blue = np.linspace(th-np.pi, th, curve_pts)
        x_blue = x + r*np.sin(th_blue)
        y_blue = y + r + r*np.cos(th_blue)
        
        # plot
        ax.fill(x_grn, y_grn, color='g', edgecolor=None, alpha=alpha)
        ax.fill(x_blue, y_blue, color='b', edgecolor=None, alpha=alpha)

    @staticmethod   
    def draw_ceiling(ax, thickness=None):
        """
        Draw ceiling for bouncing ball example
        """

        if thickness is None:
            thickness = 0.01*l

        x_ceil = [-l/2., -l/2., l/2., l/2.]
        y_ceil = [d, d+thickness, d+thickness, d]

        gray = [0.2, 0.2, 0.2]

        ax.fill(x_ceil, y_ceil, color=gray, edgecolor='k')

    @staticmethod
    def draw_paddle(ax, x, y, thickness=None):
        """
        Draw paddle for bouncing ball example
        """

        if thickness is None:
            thickness = 0.01*l

        x_paddle = [-l/2., -l/2., l/2., l/2.]
        y_paddle = [0., -thickness, -thickness, 0.]
        x_paddle = [xp + x for xp in x_paddle]
        y_paddle = [yp + y for yp in y_paddle]

        ax.fill(x_paddle, y_paddle, color='r', edgecolor='r')

    def get_fig_bounds(self, dw=None):
        """
        Get figure bounds for ball and paddle example
        """

        if dw is None:
            dw = 0.05*l

        x_min = -l/2.
        x_max = l/2.
        y_min = -l/2.
        y_max = d

        for x_k in self.x_traj:
            x_paddle = x_k[3]
            y_paddle = x_k[4]

            x_min = min(x_paddle - l/2., x_min)
            x_max = max(x_paddle + l/2., x_max)
            y_min = min(y_paddle, y_min)

        return x_min-dw, x_max+dw, y_min-dw, y_max+dw


# load modes and domains
D_arr = []
sys_arr = []
for i in range(1,8):
    # domain
    D_filename = 'D' + str(i) + '.json'
    D_arr.append(zono.from_json(D_filename))

    # system
    S_filename = 'S' + str(i) + '.json'
    sys_arr.append(AffineSystem.from_json(S_filename))

S = zono.from_json('S.json')
U = zono.from_json('U.json')

# make graph of function
Psi_i_arr = []

for i in range(len(D_arr)):
    sys_i = sys_arr[i]
    SU_i = D_arr[i]
    Psi_i_arr.append(make_graph_of_function(sys_i.A, sys_i.B, sys_i.c, SU_i))
    
Psi = zono.union_of_many(Psi_i_arr, preserve_sharpness=False)
Psi = zono.intersection_over_dims(Psi, S, [i for i in range(12, 22)])


# cost function
Q = sparse.diags([1., 1., 0.01, 1., 1., 1., 1., 0.01, 1., 1.])
R = sparse.diags([0.01, 0.001])

# number of time steps
N = 20

# initial condition
x0 = np.array([0., 0., np.pi, 0., 0., 0., 0., 0., 0., 0.])

# terminal set
xN = np.zeros(10)
SN = zono.Point(xN)
Psi_N = zono.intersection_over_dims(Psi, SN, [i for i in range(12, 22)])

# reference
xr = x0.copy()
xr_arr = [xr for _ in range(N)]

# array of graphs of functions
Psi_arr = [Psi for _ in range(N-1)]
Psi_arr.append(Psi_N)

# build and solve planning problem
settings = zono.OptSettings()
settings.t_max = 60.
settings.k_max_admm_fp_ph1 = 10000
settings.k_max_admm_fp_ph2 = 900000
settings.verbosity_interval = 1000
settings.single_threaded_admm_fp = True
settings.polish = False
settings.enable_rng_seed = True

settings.eps_prim = 0.01

if MODE == 'representative':
    settings.verbose = True
    settings.rng_seed = 4

    sol = zono.OptSolution()
    x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u = zono_planning_prob(x0, xr_arr, Psi_arr, S, U, Q, R, Q, N, settings, sol)
    x_traj_sim = simulate_pwa(x0, u_traj, sys_arr, D_arr)

    # make figure and animation
    plot_obj = BallPlot(x_traj, x_traj_sim)
    plot_obj.make_ball_paddle_fig(k_list=[0, 4, 8, 12, 16, 20], ncols=3)
    plot_obj.make_ball_paddle_animation( playback_speed=2.0)

elif MODE == 'multi-run':
    settings.verbose = False

    n_runs = 50

    iter_arr = []
    time_arr = []

    for seed in range(n_runs):
        print(f"Run {seed+1} / {n_runs}")

        settings.rng_seed = seed

        sol = zono.OptSolution()
        x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u = zono_planning_prob(x0, xr_arr, Psi_arr, S, U, Q, R, Q, N, settings, sol)

        if sol.converged:
            iter_arr.append(sol.iter)
            time_arr.append(sol.run_time)

    print(f'Number converged: {len(iter_arr)} / {n_runs}')
    print(f'Iterations: median = {np.median(iter_arr)}, min = {np.min(iter_arr)}, max = {np.max(iter_arr)}')
    print(f'Run time: median = {np.median(time_arr)}, min = {np.min(time_arr)}, max = {np.max(time_arr)}')

else:
    raise ValueError("Invalid MODE selected.")

