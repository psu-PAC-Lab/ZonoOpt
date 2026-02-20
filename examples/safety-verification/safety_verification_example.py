import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import zonoopt as zono
import subprocess

try:
    import control
except ImportError:
    print('control not found')
    exit()

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

### System model ###

# time step
dt = 0.5

# 2D double integrator dynamics
# x = [x, y, xdot, y_dot]
# u = [x_ddot, y_ddot]
A = sp.csc_matrix([[1., 0., dt, 0.],
              [0., 1., 0., dt],
              [0., 0., 1., 0.],
              [0., 0., 0., 1.]])
B = sp.csc_matrix([[0.5*dt**2, 0.],
                [0., 0.5*dt**2],
                [dt, 0.],
                [0., dt]])
nx = 4
nu = 2

# feasible state set
n_sides_circle_approx = 6
v_max = 1.0 # max velocity
Xp = zono.make_regular_zono_2D(500.0, n_sides_circle_approx)
Xv = zono.make_regular_zono_2D(v_max, n_sides_circle_approx)
S = zono.cartesian_product(Xp, Xv)

### Controller ###

# number of time steps
n_time = 24

# desired trajectory
v_traj = 0.5
phi_traj = np.linspace(np.pi/4., -3*np.pi/4., 40)
x_traj = [0.]
y_traj = [0.]
for i in range(n_time):
    x_traj.append(x_traj[-1] + v_traj*np.cos(phi_traj[i])*dt)
    y_traj.append(y_traj[-1] + v_traj*np.sin(phi_traj[i])*dt)

# LQR cost function
Q = np.diag([1., 1., 0., 0.])
R = 0.1*np.eye(2)

# LQR gain
K, slqr, elqr = control.dlqr(A.toarray(), B.toarray(), Q, R)
K = sp.csc_matrix(K)

# closed loop dynamics
Acl = A - B.dot(K)

### Initial state set ###
xc = np.array([1., 0., 0., 0.])
Xp = zono.make_regular_zono_2D(0.5, n_sides_circle_approx)
Xv = zono.make_regular_zono_2D(0.5, n_sides_circle_approx)
X0 = zono.affine_map(zono.cartesian_product(Xp, Xv), sp.eye(4), xc)

### Disturbance set ###
wp_bound = 0.01 # position noise
wv_bound = 0.2 # velocity noise
wvx0 = 0.0 # velocity bias
wvy0 = 0.5 # velocity bias
Wp = zono.make_regular_zono_2D(wp_bound, n_sides_circle_approx)
Wv = zono.make_regular_zono_2D(wv_bound, n_sides_circle_approx)
Wv = zono.affine_map(Wv, sp.eye(2), np.array([wvx0, wvy0]))
W = zono.cartesian_product(Wp, Wv)  

### unsafe set ###
V = np.array([[3.94, 0.75],
              [3.65, 1.0],
              [3.45, 1.6],
              [4.2, 1.7],
              [4.5, 0.6],
              [4.63, 1.3],
              [3.81, 1.74]])
V += np.array([0.18, 0.18])

O = zono.vrep_2_conzono(V)

### Get reachable sets ###
X = X0
X_arr = [X]
for k in range(n_time):

    xr = np.array([x_traj[k], y_traj[k], 0., 0.])
    uff = K.dot(xr)

    XWS = zono.cartesian_product(zono.cartesian_product(X, W), S)
    AImI = sp.hstack((Acl, sp.eye(nx), -sp.eye(nx)))
    XWS_int = zono.intersection(XWS, zono.Point(-B.dot(uff)), AImI)
    X = zono.project_onto_dims(XWS_int, [i for i in range(2*nx, 3*nx)])

    X_arr.append(X)

### Check for safety of reachable sets ###
sol = zono.OptSolution()
settings = zono.OptSettings()
settings.k_inf_check = 1
settings.rho = 1.

iter_arr = []
soltime_arr = []
for X in X_arr:
    if zono.intersection_over_dims(X, O, [0,1]).is_empty(solution=sol, settings=settings):
        print(f'Safe: iter = {sol.iter}, sol time = {sol.run_time}')
    else:
        print(f'Unsafe: iter = {sol.iter}, sol time = {sol.run_time}')
    iter_arr.append(sol.iter)
    soltime_arr.append(sol.run_time)
iter_arr = np.array(iter_arr)
soltime_arr = np.array(soltime_arr)

print(f'Max solution time [ms]: {1000.*np.max(soltime_arr)}')
print(f'Min solution time [ms]: {1000.*np.min(soltime_arr)}')
print(f'Average solution time [ms]: {1000.*np.mean(soltime_arr)}')

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
figsize = (245.71 * inches_per_pt, 0.8*245.71 * inches_per_pt)  # Convert pt to inches


with plt.rc_context(rc_context):

    fig = plt.figure(constrained_layout=True, figsize=figsize)
    ax = fig.add_subplot(111)

    # colorbars
    cmap = plt.get_cmap('cool', np.max(iter_arr)-np.min(iter_arr)+1)
    norm = mpl.colors.Normalize(vmin=np.min(iter_arr), vmax=np.max(iter_arr))
    sm = plt.cm.ScalarMappable(cmap='cool', norm=norm)
    sm.set_array([])
    
    xr_obj = ax.plot(x_traj, y_traj, '.r')
    for k, X in enumerate(X_arr):
        color = cmap(norm(iter_arr[k]))
        if k == 0:
            X_obj = zono.plot(zono.project_onto_dims(X, [0,1]), ax, color=color, edgecolor='k', label='Reachable set')
        else:
            zono.plot(zono.project_onto_dims(X, [0,1]), ax, color=color, edgecolor='k')
    O_obj = zono.plot(O, ax, color='g', edgecolor='k')
    ax.axis('equal')
    plt.colorbar(sm, ax=ax, label='ADMM iterations')
    ax.grid(alpha=0.2)
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$y$ [m]')

    ax.legend(handles=[X_obj[0], O_obj[0], xr_obj[0]], 
              labels=[r'$\mathcal{X}_k$', r'$\mathcal{O}$', r'$\mathbf{x}_k^r$'])

    if is_latex_installed():
        plt.savefig('reachability_traj.pgf')

    plt.show()