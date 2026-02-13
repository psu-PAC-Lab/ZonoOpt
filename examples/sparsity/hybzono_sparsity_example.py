import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import zonoopt as zono
import subprocess

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def plot_reach_sets(X_arr, title=None):

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
    figsize = (245.71 * inches_per_pt, 0.8*245.71 * inches_per_pt)  # Convert pt to inches

    with plt.rc_context(rc_context):

        # create figure
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        ax = fig.add_subplot(1,1,1)

        # colorbars
        cmap = plt.get_cmap('cool', len(X_arr))
        norm = mpl.colors.Normalize(vmin=0, vmax=len(X_arr)-1)
        sm = plt.cm.ScalarMappable(cmap='cool', norm=norm)
        sm.set_array([])

        # plot sets
        for k, X in enumerate(X_arr):
            print(f'Plotting reachability step {k}')
            color = cmap(norm(k))
            zono.plot(X, ax=ax, color=color)

        ax.axis('equal')
        plt.colorbar(sm, ax=ax, label=r'$k$', ticks=range(len(X_arr)))
        if title is not None:
            ax.set_title(title)
        ax.grid(alpha=0.2)
        ax.set_xlabel(r'$x_1$')
        ax.set_ylabel(r'$x_2$')

        if is_latex_installed():
            plt.savefig('two_equilibrium_reach_sets.pgf')

        plt.show()

def plot_convex_relaxations(X_arr_mld, X_arr_gof_cond, X_arr_gof_sharp):

    k = 5  # time step to plot

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
    figsize = (245.71 * inches_per_pt, 0.8*245.71 * inches_per_pt)  # Convert pt to inches

    with plt.rc_context(rc_context):

        # create figure
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        ax = fig.add_subplot(1,1,1)
    
        # plot convex relaxations
        soft_frac = 0.3
        soft_red = (1., soft_frac, soft_frac)
        soft_blue = (soft_frac, soft_frac, 1.)
        soft_green = (soft_frac, 1., soft_frac)

        h = []
        h.append(zono.plot(X_arr_mld[k].convex_relaxation(), ax=ax, color=soft_red, alpha=0.3, edgecolor='k')[0])
        h.append(zono.plot(X_arr_gof_cond[k].convex_relaxation(), ax=ax, color=soft_blue, alpha=0.3, edgecolor='k')[0])
        h.append(zono.plot(X_arr_gof_sharp[k].convex_relaxation(), ax=ax, color=soft_green, alpha=0.3, edgecolor='k')[0])
        h.append(zono.plot(X_arr_gof_cond[k], ax=ax, color='k', alpha=1.)[0])
        
        ax.set_aspect('equal', adjustable='box')
        ax.legend(handles=h,
                    labels=(r'$\mathit{CR}(\mathcal{X}_5)$ MLD', 
                          r'$\mathit{CR}(\mathcal{X}_5)$ PWA -- Cond. Union', 
                          r'$\mathit{CR}(\mathcal{X}_5)$ PWA -- Sharp Union', 
                          r'$\mathcal{X}_5$'), fontsize=textwidth_pt)
        ax.grid(alpha=0.2)
        ax.set_xlabel(r'$x_1$')
        ax.set_ylabel(r'$x_2$')

        # ax.set_ylim((-1.5, 4.5))

        if is_latex_installed():
            plt.savefig('two_equilibrium_convex_relaxations.pgf')

        plt.show()

def print_complexity(X):
    # [nGc, nGb, nC, G nnz, A nnz]
    n = X.get_n()
    nGc = X.get_nGc()
    nGb = X.get_nGb()
    nC = X.get_nC()
    G = X.get_G()
    A = X.get_A()
    return f'[{n} {nGc} {nGb} {nC} {G.nnz} {A.nnz}]'


# two equilibrium example from Bird, Trevor J., et al. "Hybrid zonotopes: A new set representation for reachability analysis of mixed logical dynamical systems." Automatica 154 (2023): 111107.
n = 2
nu = 0
N = 15

Ad1 = np.array([[0.75, 0.25],
                [-0.25, 0.75]])
Bd1 = np.zeros((2,0))
fd1 = np.array([-0.25, -0.25])
Ad2 = np.array([[0.75, -0.25],
                [0.25, 0.75]])
Bd2 = np.zeros((2,0))
fd2 = np.array([0.25, -0.25])

# feasible sets for states / inputs
SU1 = zono.interval_2_zono(zono.Box([-2., -1.], [0., 3]))
SU2 = zono.interval_2_zono(zono.Box([0., -1.], [2., 3.]))
S = zono.interval_2_zono(zono.Box([-2., -1.], [2., 3.]))

# input set
U = zono.EmptySet(0)

# MLD formulation (ref https://github.com/ESCL-at-UTD/zonoLAB/blob/main/examples/MLD_Reachability/exampleTwoEquilibrium.m, generated by Hysdel)
ne = 10
A_mld = np.zeros((n,n))
B_u_mld = np.zeros((n,nu))
B_w_mld = np.array([[1., 0., 0.],
                    [0., 1., 0.]])
B_aff_mld = np.zeros(n)
E_x_mld = np.array([
    [1., 0.],
    [-1., 0.],
    [0.75, -0.25],
    [-0.75, 0.25],
    [0.75, 0.25],
    [-0.75, -0.25],
    [0.25, 0.75],
    [-0.25, -0.75],
    [-0.25, 0.75],
    [0.25, -0.75]
])
E_u_mld = np.zeros((ne,nu))
E_w_mld = np.array([
    [0., 0., -5.],
    [0., 0., 5.],
    [-1., 0., 10.5],
    [1., 0., 9.5],
    [-1., 0., -9.5],
    [1., 0., -10.5],
    [0., -1., 10.], 
    [0., 1., 10.],
    [0., -1., -10.],
    [0., 1., -10.]
])
E_aff_mld = np.array([0., 5., 10.25, 9.75, 0.25, -0.25, 10.25, 9.75, 0.25, -0.25])

AE_x = np.vstack((A_mld, E_x_mld))
BE_u = np.vstack((B_u_mld, E_u_mld))
BE_w = np.vstack((B_w_mld, E_w_mld))
BO_aff = np.hstack((B_aff_mld, np.zeros(ne)))
OI_aff = np.hstack((np.zeros((ne,n)), np.identity(ne)))


# initial set
X0 = zono.Zono(np.array([[0.25, -0.19],
                         [0.19, 0.25]]), 
                np.array([-1.31, 2.55]))

# auxiliary variable set
W_cvx = zono.interval_2_zono(zono.Box([-5.25, -5.25, 0.], [5.25, 4.74, 1.]))
W = zono.HybZono(W_cvx.get_G()[:, 0:2], W_cvx.get_G()[:,2], W_cvx.get_c(),
                 np.zeros((0,2)), np.zeros((0,1)), np.zeros(0))                           


### standard, unconstrained reachability

# MLD reachability
X = X0.copy()
X_arr_mld = [X]
V = zono.minkowski_sum(zono.minkowski_sum(zono.affine_map(U, BE_u), zono.affine_map(W, BE_w)),
                        zono.Point(BO_aff))
for _ in range(N):
    X = zono.minkowski_sum(zono.affine_map(X, AE_x), V)
    X = zono.halfspace_intersection(X, np.identity(ne), E_aff_mld, OI_aff)
    X = zono.project_onto_dims(X, [i for i in range(n)])
    X_arr_mld.append(X)

# plot_reach_sets(X_arr_mld, 'MLD Reachable Sets')

# GOFs for each mode
IIAB_1 = np.vstack( (np.eye(n+nu), np.hstack( (Ad1, Bd1) )))
fd1_ext = np.hstack( (np.zeros(n+nu), fd1) )
Psi_1 = zono.minkowski_sum(zono.affine_map(SU1, IIAB_1), zono.Point(fd1_ext))
IIAB_2 = np.vstack( (np.eye(n+nu), np.hstack( (Ad2, Bd2) )))
fd2_ext = np.hstack( (np.zeros(n+nu), fd2) )
Psi_2 = zono.minkowski_sum(zono.affine_map(SU2, IIAB_2), zono.Point(fd2_ext))

# GOF reachability - condensed union
Psi = zono.union_of_many((Psi_1, Psi_2), preserve_sharpness=False) 
X = X0.copy()
X_arr_gof_cond = [X]
for _ in range(N):
    X = zono.project_onto_dims(zono.intersection_over_dims(Psi, X, [i for i in range(n)]), [i for i in range(n+nu, n+nu+n)] )
    X_arr_gof_cond.append(X)

plot_reach_sets(X_arr_gof_cond)

# GOF reachability - sharp union
Psi = zono.union_of_many((Psi_1, Psi_2), preserve_sharpness=True)
X = X0.copy()
X_arr_gof_sharp = [X]
for _ in range(N):
    X = zono.project_onto_dims(zono.intersection_over_dims(Psi, X, [i for i in range(n)]), [i for i in range(n+nu, n+nu+n)] )
    X_arr_gof_sharp.append(X)

# GOF reachability - zonotope union
Psi = zono.zono_union_2_hybzono((Psi_1, Psi_2))
X = X0.copy()
X_arr_gof_zono = [X]
for _ in range(N):
    X = zono.project_onto_dims(zono.intersection_over_dims(Psi, X, [i for i in range(n)]), [i for i in range(n+nu, n+nu+n)] )
    X_arr_gof_zono.append(X)

# plot convex relaxations
plot_convex_relaxations(X_arr_mld, X_arr_gof_cond, X_arr_gof_sharp)


### lifted, constrained reachability

# GOF reachability - condensed union
Psi = zono.union_of_many((Psi_1, Psi_2), preserve_sharpness=False)
Psi = zono.intersection_over_dims(Psi, S, [i for i in range(n+nu, n+nu+n)])
Z = X0.copy()
X_arr_gof_cond_lifted = [Z]
for _ in range(N):
    ZUS = zono.cartesian_product(Z, S) # U is empty so leaving out
    Z = zono.intersection_over_dims(ZUS, Psi, [i for i in range(ZUS.get_n()-n-nu-n, ZUS.get_n())])
    X_arr_gof_cond_lifted.append(zono.project_onto_dims(Z, [i for i in range(Z.get_n()-n, Z.get_n())]))
Z_gof_cond = Z.copy()

# plot_reach_sets(X_arr_gof_cond_lifted)

# GOF reachability - sharp union
Psi = zono.union_of_many((Psi_1, Psi_2), preserve_sharpness=True)
Psi = zono.intersection_over_dims(Psi, S, [i for i in range(n+nu, n+nu+n)])
Z = X0.copy()
X_arr_gof_sharp_lifted = [Z]
for _ in range(N):
    ZUS = zono.cartesian_product(Z, S) # U is empty so leaving out
    Z = zono.intersection_over_dims(ZUS, Psi, [i for i in range(ZUS.get_n()-n-nu-n, ZUS.get_n())])
    X_arr_gof_sharp_lifted.append(zono.project_onto_dims(Z, [i for i in range(Z.get_n()-n, Z.get_n())]))
Z_gof_sharp = Z.copy()

# GOF reachability - zonotope union
Psi = zono.zono_union_2_hybzono((Psi_1, Psi_2))
Psi = zono.intersection_over_dims(Psi, S, [i for i in range(n+nu, n+nu+n)])
Z = X0.copy()
X_arr_gof_zono_lifted = [Z]
for _ in range(N):
    ZUS = zono.cartesian_product(Z, S) # U is empty so leaving out
    Z = zono.intersection_over_dims(ZUS, Psi, [i for i in range(ZUS.get_n()-n-nu-n, ZUS.get_n())])
    X_arr_gof_zono_lifted.append(zono.project_onto_dims(Z, [i for i in range(Z.get_n()-n, Z.get_n())]))
Z_gof_zono = Z.copy()



# table for sparsity / complexity
# [nGc, nGb, nC, Gc sparsity, Gb sparsity, Ac sparsity, Ab sparsity]
print('                           [n, nGc, nGb, nC, G nnz, A nnz]')
print('MLD:                       ' + print_complexity(X_arr_mld[-1]))
print('GOF -- condensed union:    ' + print_complexity(X_arr_gof_cond[-1]))
print('GOF -- sharp union:        ' + print_complexity(X_arr_gof_sharp[-1]))
# print('GOF -- zonotope union:      ' + print_complexity(X_arr_gof_zono[-1]))
print('GOF lifted + constrained -- condensed:    ' + print_complexity(Z_gof_cond))
print('GOF lifted + constrained -- sharp:        ' + print_complexity(Z_gof_sharp))

    
# # plot sparsity
# fig = plt.figure(constrained_layout=True)

# # MLD
# ax = fig.add_subplot(2, 3, 1)
# G = X_arr_mld[-1].get_G()
# ax.spy(G, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'MLD G: nnz = {G.nnz}')

# ax = fig.add_subplot(2, 3, 4)
# A = X_arr_mld[-1].get_A()
# ax.spy(A, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'MLD A: nnz = {A.nnz}')

# # GOF
# ax = fig.add_subplot(2, 3, 2)
# G = X_arr_gof_cond[-1].get_G()
# ax.spy(G, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'GOF cond. G: nnz = {G.nnz}')

# ax = fig.add_subplot(2, 3, 5)
# A = X_arr_gof_cond[-1].get_A()
# ax.spy(A, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'GOF cond. A: nnz = {A.nnz}')

# # original GOF
# ax = fig.add_subplot(2, 3, 3)
# G = X_arr_gof_sharp[-1].get_G()
# ax.spy(G, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'GOF sharp G: nnz = {G.nnz}')

# ax = fig.add_subplot(2, 3, 6)
# A = X_arr_gof_sharp[-1].get_A()
# ax.spy(A, color=(0.5,0.5,0.5), aspect='equal', markersize=3)
# ax.set_title(f'GOF sharp A: nnz = {A.nnz}')

# plt.show()