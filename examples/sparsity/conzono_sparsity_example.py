import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import gridspec
import zonoopt as zono
import subprocess

# second order system
# x = [z, zdot]
# zddot + 2*xi*wn*zdot + wn^2*z = u
dt = 0.1
wn = 0.3
xi = 0.7
A = sp.csc_matrix([[1., dt], 
                   [-wn**2*dt, 1-2*xi*wn*dt]])
B = sp.csc_matrix([[0.], 
                   [dt]])

nx = 2
nu = 1

# initial set
X0 = zono.Zono(sp.diags([0.01, 0.01]), np.array([0., 0.5]))

# input set
U = zono.Zono(sp.diags([1.0]), np.array([0.]))

# state constraints
S = zono.Zono(sp.diags([1.0, 1.0]), np.array([0., 0.]))

# number of time steps
N = 15

### Reachability analysis with linear map, minkowski sum, and intersection ###
X1 = X0
X1_arr = [X1]
for k in range(N):
    X1 = zono.affine_map(X1, A)
    X1 = zono.minkowski_sum(X1, zono.affine_map(U, B))
    X1 = zono.intersection(X1, S)
    X1_arr.append(X1)

### Reachability analysis with state update set ###
IAB = sp.vstack(( sp.eye(nx+nu), sp.hstack((A, B)) ))
Psi = zono.affine_map(zono.cartesian_product(S, U), IAB)
Psibar = zono.intersection_over_dims(Psi, S, [i for i in range(nx+nu, nx+nu+nx)])

X2 = X0
X2_arr = [X2]
for k in range(N):
    XU = zono.cartesian_product(X2, U)
    PsiRU = zono.intersection_over_dims(Psibar, XU, [i for i in range(nx+nu)])
    X2 = zono.project_onto_dims(PsiRU, [i for i in range(nx+nu, nx+nu+nx)])
    X2_arr.append(X2)

### Reachability analysis with method from paper ###
X3 = X0
X3_arr = [X3]
for k in range(N):
    XUS = zono.cartesian_product(zono.cartesian_product(X3, U), S)
    ABmI = sp.hstack((A, B, -sp.eye(nx)))
    XUS_int = zono.intersection(XUS, zono.Point(np.zeros(nx)), ABmI)
    X3 = zono.project_onto_dims(XUS_int, [i for i in range(nx+nu, nx+nu+nx)])
    X3_arr.append(X3)


### plot
def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

# print sparsity
print(f'X1 sparsity: nnz(G) = {X1.get_G().nnz}, density(G) = {X1.get_G().nnz/(X1.get_n()*X1.get_nG())}, nnz(A) = {X1.get_A().nnz}, density(A) = {X1.get_A().nnz/(X1.get_nC()*X1.get_nG())}')
print(f'X2 sparsity: nnz(G) = {X2.get_G().nnz}, density(G) = {X2.get_G().nnz/(X2.get_n()*X2.get_nG())}, nnz(A) = {X2.get_A().nnz}, density(A) = {X2.get_A().nnz/(X2.get_nC()*X2.get_nG())}')
print(f'X3 sparsity: nnz(G) = {X3.get_G().nnz}, density(G) = {X3.get_G().nnz/(X3.get_n()*X3.get_nG())}, nnz(A) = {X3.get_A().nnz}, density(A) = {X3.get_A().nnz/(X3.get_nC()*X3.get_nG())}')


# set up LaTeX rendering
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

# make plot
X_arr_arr = [X1_arr.copy(), X2_arr.copy(), X3_arr.copy()]
col_titles = [r'Eq.~\ref{eq:reachability-trad}' + '\n', 
              r'Eq.~\ref{eq:reachability-state-update}' + '\n' + r'$G$ sparsity', 
              r'Eq.~\ref{eq:reachability-alt-int}' + '\n']


# figure size
figwidth_pt = 505.89
inches_per_pt = 1 / 72.27
figsize = (figwidth_pt * inches_per_pt, 0.45*figwidth_pt * inches_per_pt)  # Convert pt to inches

with plt.rc_context(rc_context):
    fig = plt.figure(constrained_layout=True, figsize=figsize)

    width_col_pt = (figwidth_pt / 3.) - textwidth_pt
    gs = fig.add_gridspec(nrows=2, ncols=3, figure=fig, width_ratios=[1,1,1], height_ratios=[1,3])

    for row in range(2):
        for col in range(3):

            # sets
            X_arr = X_arr_arr[col]
            XN = X_arr_arr[col][N]

            ax = fig.add_subplot(gs[row,col])

            if row == 0:

                # G sparsity
                ax.spy(XN.get_G(), color=(0.5,0.5,0.5), aspect='equal', markersize=3)
                ax.set_yticks(ticks=[0,1],labels=['', '1'])
                
                ax.set_title(col_titles[col], fontsize=textwidth_pt)

            else:

                # A sparsity
                ax.spy(XN.get_A(), color=(0.5,0.5,0.5), aspect='equal', markersize=3)
                
                if col == 1:
                    ax.set_title(r'$A$ sparsity', fontsize=textwidth_pt)

    # axis labels
    fig.supxlabel(r'columns', fontsize=textwidth_pt)
    fig.supylabel(r'rows', fontsize=textwidth_pt)

    # save
    if is_latex_installed():
        plt.savefig('reachability_sparsity.pgf')

    plt.show()


# reachability plot
figwidth_pt = 245.71
figsize = (figwidth_pt * inches_per_pt, 0.7*figwidth_pt * inches_per_pt)  # Convert pt to inches

with plt.rc_context(rc_context):
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    ax = fig.add_subplot(111)

    for i in range(N):
        zono.plot(X_arr[i], ax, color='b', alpha=0.1)
    zono.plot(XN, ax, color='k', alpha=0.2)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v$')

    if is_latex_installed():
        plt.savefig('reachability_reach_sets.pgf')

    plt.show()
       