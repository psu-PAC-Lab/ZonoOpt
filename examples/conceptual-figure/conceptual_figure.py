import zonoopt as zono
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import subprocess

# 2D second order system
dt = 1.0
wn = 0.325
xi = 0.7
A = sparse.csc_matrix([[1., 0., dt, 0.], 
                     [0., 1., 0., dt],
                   [-wn**2*dt, 0., 1-2*xi*wn*dt, 0.],
                   [0., -wn**2*dt, 0., 1-2*xi*wn*dt]])
B = sparse.csc_matrix([[0., 0.],
                    [0., 0.], 
                    [dt, 0.], 
                    [0., dt]])
nx = 4
nu = 2

# state feasible set
S_pos = zono.interval_2_zono(zono.Box([0., 0.], [10., 10.]))
S_vel = zono.minkowski_sum(
    zono.affine_map(zono.make_regular_zono_2D(0.2, 8), sparse.csc_matrix([[1., 0.], [0., 3.]])),
    zono.Point(np.array([1., 0.]))
)
S = zono.cartesian_product(S_pos, S_vel)

# input feasible set
U = zono.make_regular_zono_2D(1., 8)

# initial set
X0_pos = zono.make_regular_zono_2D(0.1, 8, c=np.array([0.5, 2.]))
X0_vel = zono.make_regular_zono_2D(0.1, 8, c=np.array([1.0, -0.3]))
X0 = zono.cartesian_product(X0_pos, X0_vel)

# reference (constant)
xr = np.array([10., 5., 0., 0.])

# cost matrices
Q = sparse.diags([0.1, 0.1, 0.0, 0.0])
R = sparse.diags([5., 5.])

### begin reachability analysis ###

N = 7  # number of time steps
idx = 0
idx_x = []
idx_u = []

idx_x.append([j for j in range(idx, idx+nx)])
idx += nx

# init cost
P = Q
q = np.zeros(nx)

# init set
Z = X0

for k in range(N):

    # control
    Z = zono.cartesian_product(Z, U)

    idx_u.append([j for j in range(idx, idx+nu)])
    idx += nu
    
    # dynamics
    nZ = Z.get_n()
    Z = zono.cartesian_product(Z, S)
    ABmI = sparse.hstack((sparse.csc_matrix((nx, nZ-(nx+nu))), A, B, -sparse.eye(nx)))
    Z = zono.intersection(Z, zono.Point(np.zeros(nx)), ABmI)

    idx_x.append([j for j in range(idx, idx+nx)])
    idx += nx

    # cost
    P = sparse.block_diag((P, R, Q))
    q = np.hstack([q, np.zeros(nu), -Q.dot(xr)])

### solve optimization problem ###
z_sol = Z.optimize_over(P, q)

### plot ###
def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

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


with plt.rc_context(rc_context):

    fig = plt.figure(constrained_layout=True, figsize=figsize)
    ax = fig.add_subplot(111)

    for i in range(N):
        X = zono.project_onto_dims(Z, idx_x[i][0:2])
        zono.plot(X, alpha=0.2, color='blue')
        ax.plot(z_sol[idx_x[i][0]], z_sol[idx_x[i][1]], 'or', markersize=4)

    ax.grid(alpha=0.2)
    ax.legend(['Reachable sets', 'Optimal trajectory'], loc='upper left')

    ax.arrow(0.75, 0.5, 1.5, 0.0, width=0.03, color='m')

    if is_latex_installed():
        plt.savefig('zono_opt_conceptual.pgf')

    plt.show()