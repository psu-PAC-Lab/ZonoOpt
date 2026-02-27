import zonoopt as zono
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import subprocess

# PWA system
nx = 4
nu = 2

dt = 1.0
wn = 0.3
xi = 0.3

A1 = np.array([[1., 0., dt, 0.], 
                     [0., 1., 0., dt],
                   [-wn**2*dt, 0., 1-2*xi*wn*dt, 0.],
                   [0., -wn**2*dt, 0., 1-2*xi*wn*dt]])
B1 = np.array([[0., 0.],
                    [0., 0.], 
                    [wn**2*dt, 0.], 
                    [0., wn**2*dt]])
f1 = np.array([0.5, 0., 0., -0.1])

A2 = A1
B2 = B1
f2 = np.array([0.5, 0., 0., 0*0.05])

# domains
Sp1 = zono.interval_2_zono(zono.Box([-5., -5.], [15., 0.]))
Sp2 = zono.interval_2_zono(zono.Box([-5., 0.], [15., 5.]))

Sv = zono.make_regular_zono_2D(1., 6)

U1 = zono.vrep_2_conzono(np.array([[-0.08, 0.15],
                                 [-0.15, 0.1],
                                 [0.12, -0.16],
                                 [0.15, -0.1]]))
U2 = zono.vrep_2_conzono(np.array([[-0.2, 0.2],
                                    [-0.23, 0.17],
                                    [0.24, -0.2],
                                    [0.16, -0.15]]))

SU1 = zono.cartesian_product(zono.cartesian_product(Sp1, Sv), U1)
SU2 = zono.cartesian_product(zono.cartesian_product(Sp2, Sv), U2)

# overall bounds
Sbar = zono.interval_2_zono(zono.Box([-5., -5., -1., -1.], [15., 5., 1., 1.]))
U = zono.union_of_many((U1, U2))
Ubar = zono.interval_2_zono(U.bounding_box())

# GOFs
IIAB_1 = np.vstack( (np.eye(nx+nu), np.hstack( (A1, B1) )))
fd1_ext = np.hstack( (np.zeros(nx+nu), f1) )
Psi_1 = zono.minkowski_sum(zono.affine_map(SU1, IIAB_1), zono.Point(fd1_ext))
IIAB_2 = np.vstack( (np.eye(nx+nu), np.hstack( (A2, B2) )))
fd2_ext = np.hstack( (np.zeros(nx+nu), f2) )
Psi_2 = zono.minkowski_sum(zono.affine_map(SU2, IIAB_2), zono.Point(fd2_ext))

Psi = zono.union_of_many((Psi_1, Psi_2))

# initial set
X0 = zono.cartesian_product(zono.make_regular_zono_2D(0.1, 8, c=[-2., -0.5]),
                            zono.make_regular_zono_2D(0.01, 4, c=[-0.5, 0.3]))

# reference (constant)
xr = np.array([10., -0.15, 0., 0.])

# cost matrices
Q = sparse.diags([0.1, 1., 0., 0.])
QN = 7*Q
R = 0.01*sparse.diags([0.1, 0.1])

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
c = 0.

# init set
Z = X0

for k in range(N):

    # lift
    Z = zono.cartesian_product(Z, Ubar)
    Z = zono.cartesian_product(Z, Sbar)

    idx_u.append([j for j in range(idx, idx+nu)])
    idx += nu

    idx_x.append([j for j in range(idx, idx+nx)])
    idx += nx
    
    # dynamics and constraints
    nZ = Z.get_n()
    Z = zono.intersection_over_dims(Z, Psi, [i for i in range(nZ - Psi.get_n(), nZ)])

    # cost
    if k == N-1:
        P = sparse.block_diag((P, R, QN))
        q = np.hstack([q, np.zeros(nu), -QN.dot(xr)])
        c += 0.5*xr.dot(QN.dot(xr))
    else:
        P = sparse.block_diag((P, R, Q))
        q = np.hstack([q, np.zeros(nu), -Q.dot(xr)])
        c += 0.5*xr.dot(Q.dot(xr))

### solve optimization problem ###
settings = zono.OptSettings()
settings.single_threaded_admm_fp = True
settings.rng_seed = 0
settings.enable_rng_seed = True
z_sol = Z.optimize_over(P, q, settings=settings)

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

    h = [None, None]
    for i in range(N+1):
        X = zono.project_onto_dims(Z, idx_x[i][0:2])
        h[0] = zono.plot(X, color=(205./255., 201./255., 255./255.), edgecolor=None)[0]
        h[1] = ax.plot(z_sol[idx_x[i][0]], z_sol[idx_x[i][1]], 'or', markersize=4)[0]

    ax.grid(alpha=0.2)
    ax.legend(h, ['Reachable sets', 'Trajectory plan'], loc='upper left')

    ax.arrow(-1.25, -0.5, 1.5, 0.0, width=0.02, color='m')
    ax.set_ylim((-0.9, 0.7))

    if is_latex_installed():
        plt.savefig('zono_opt_conceptual_hybrid.pgf')

    plt.show()