import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import zonoopt as zono
import subprocess

# when approximating a circle as a regular zonotope, use this many sides
n_sides_circle_approx = 6

def mhe(X_Nm1, W, V, A, B, u_vec, y_vec, Qinv, Rinv, S, settings=zono.OptSettings()):
    """Solves constrained zonotope representation of the MHE optimization problem.
    X_Nm1: set of all possible states at time N-1.
    W: disturbance set.
    V: measurement noise set.
    A: state transition matrix.
    B: control input matrix.
    u_vec: control inputs.
    y_vec: output measurements.
    Qinv: process noise covariance.
    Rinv: measurement noise covariance.
    S: feasible state set.
    """

    # dimensions
    nx = A.shape[1]

    # number of time steps
    N = len(y_vec)

    # init indexing
    idx_x = []
    idx_w = []
    idx = 0

    # init conzono
    Z = X_Nm1

    # indexing
    idx_x.append([j for j in range(idx, nx)])
    idx += nx

    # init cost
    P = sp.csc_matrix((nx, nx))
    q = np.zeros(nx)

    # loop through time steps to build up conzono and cost
    for k in range(-N, 0):

        # control input
        u = zono.Point(u_vec[k])
        mBu = zono.affine_map(u, -B)

        # feasible set of states from S and measurement
        mV = zono.affine_map(V, -sp.eye(V.get_n()))
        y_mV = zono.minkowski_sum(zono.Point(y_vec[k]), mV)
        Sy = zono.intersection(S, y_mV)

        # propagate state and enforce constraints
        Z = zono.cartesian_product(zono.cartesian_product(Z, W), Sy)
        AImI = sp.hstack((sp.csc_matrix((nx, Z.get_n()-3*nx)), A, sp.eye(nx), -sp.eye(nx)))
        Z = zono.intersection(Z, mBu, AImI)

        # cost
        P = sp.block_diag((P, Qinv, Rinv))
        q = np.hstack((q, np.zeros(nx), -Rinv.dot(y_vec[k]).flatten()))

        # indexing
        idx_w.append([j for j in range(idx, idx+nx)])
        idx_x.append([j for j in range(idx+nx, idx+2*nx)])
        idx += 2*nx

    # set of all possible states at time k=0
    X0 = zono.project_onto_dims(Z, idx_x[N])

    # solve
    sol = zono.OptSolution() # pass by reference
    xopt = Z.optimize_over(P, q, solution=sol, settings=settings)

    # check that solution is feasible
    if sol.infeasible:
        raise ValueError("MHE problem is infeasible.")
    
    # get state estimate for current time step
    xhat = xopt[idx_x[N]]

    return xhat, sol, X0
    
def state_update(X, A, B, W, V, u, y, S):
    """Get possible states at k+1.
    X: set of all possible states at time k.
    A: state transition matrix.
    B: control input matrix.
    W: disturbance set.
    V: measurement noise set.
    u: control input.
    y: output measurement.
    S: feasible state set.
    """

    nx = A.shape[1]

    # control input
    mBu = zono.affine_map(zono.Point(u), -B)

    # feasible set of states from S and measurement
    mV = zono.affine_map(V, -sp.eye(V.get_n()))
    y_mV = zono.minkowski_sum(zono.Point(y), mV)
    Sy = zono.intersection(S, y_mV)

    # propagate state and enforce constraints
    Z = zono.cartesian_product(zono.cartesian_product(X, W), Sy)
    AImI = sp.hstack((A, sp.eye(nx), -sp.eye(nx)))
    X_kp1 = zono.project_onto_dims(zono.intersection(Z, mBu, AImI), [i for i in range(2*nx, 3*nx)])

    return X_kp1

def get_noise(sigma_p, sigma_v):
    """Get noise.
    sigma_p: position noise.
    sigma_v: velocity noise.
    """

    p_mag = np.random.normal(0, sigma_p)
    if np.abs(p_mag) > 2*sigma_p:
        p_mag = 2*sigma_p*np.sign(p_mag)
    p_dir = np.random.uniform(0, 2*np.pi)
    p = np.array([p_mag*np.cos(p_dir), p_mag*np.sin(p_dir)])
    v_mag = np.random.normal(0, sigma_v)
    if np.abs(v_mag) > 2*sigma_v:
        v_mag = 2*sigma_v*np.sign(v_mag)
    v_dir = np.random.uniform(0, 2*np.pi)
    v = np.array([v_mag*np.cos(v_dir), v_mag*np.sin(v_dir)])
    
    return np.hstack((p, v))

def is_latex_installed():
    try:
        subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


### System model ###

# time step
dt = 1.0

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

# feasible state set
v_max = 1.0 # max velocity
Xp = zono.make_regular_zono_2D(500.0, n_sides_circle_approx, outer_approx=True)
Xv = zono.make_regular_zono_2D(v_max, n_sides_circle_approx, outer_approx=True)
S = zono.cartesian_product(Xp, Xv)

# process noise set
sigma_wp = 0.001 # position noise
sigma_wv = 0.01 # velocity noise
Wp = zono.make_regular_zono_2D(2*sigma_wp, n_sides_circle_approx, outer_approx=True)
Wv = zono.make_regular_zono_2D(2*sigma_wv, n_sides_circle_approx, outer_approx=True)
W = zono.cartesian_product(Wp, Wv)  

# measurement noise set
sigma_vp = 0.5 # position measurement noise
sigma_vv = 0.2 # velocity measurement noise
Vp = zono.make_regular_zono_2D(2*sigma_vp, n_sides_circle_approx, outer_approx=True)
Vv = zono.make_regular_zono_2D(2*sigma_vv, n_sides_circle_approx, outer_approx=True)
V = zono.cartesian_product(Vp, Vv)

### MHE parameters ###
Qinv = sp.diags([1/sigma_wp**2, 1/sigma_wp**2, 1/sigma_wv**2, 1/sigma_wv**2]) # process noise covariance
Rinv = sp.diags([1/sigma_vp**2, 1/sigma_vp**2, 1/sigma_vv**2, 1/sigma_vv**2]) # measurement noise covariance

### Simulate ###

# initial state set
X0_p = zono.make_regular_zono_2D(2.0, n_sides_circle_approx, outer_approx=True, c=np.array([-4., 1.]))
X0_v = zono.make_regular_zono_2D(1.0, n_sides_circle_approx, outer_approx=True)
X0 = zono.cartesian_product(X0_p, X0_v)

# true initial state
x0 = np.array([-5., 2., -0.6, 0.2])

# number of time steps
n_sim = 40

# horizon length
N = 15

# rng seed
np.random.seed(1)

# solver settings
settings = zono.OptSettings()
settings.rho = 1.
settings.inf_norm_conv = False
settings.eps_prim = 1e-2
settings.eps_dual = 1e-2

# load in u data
u_data = np.loadtxt('u_data.txt')

# possible state set at k=N-1
X_Nm1 = X0 # init

# simulate first N time steps
x = x0
y_sim = []
u_sim = []
x_sim = [x]
xhat_sim = []
sol_times = []
startup_times = []
iters = []
X0_arr = [X_Nm1]

# first time step

# control input
u = u_data[0,:]
u_sim.append(u)

# process noise
w = get_noise(sigma_wp, sigma_wv)

# state update
x = A.dot(x) + B.dot(u) + w
x_sim.append(x)

# simulate up to horizon length
for k in range(1, N):

    # measurement
    v = get_noise(sigma_vp, sigma_vv)
    y = x + v
    y_sim.append(y)

    # MHE
    u_vec = u_sim
    y_vec = y_sim

    try:
        xhat, sol, X0 = mhe(X_Nm1, W, V, A, B, u_vec, y_vec, Qinv, Rinv, S, settings=settings)
    except ValueError as e:
        print(e)
        break

    xhat_sim.append(xhat)
    sol_times.append(sol.run_time)
    startup_times.append(sol.startup_time)
    iters.append(sol.iter)
    print(f'k = {k}, sol_time = {sol.run_time}, iter = {sol.iter}, |x-xhat|_2 = {np.linalg.norm(x-xhat)}')

    # control input
    if k < u_data.shape[0]:
        u = u_data[k,:]
    else:
        u = u_data[-1,:]
    u_sim.append(u)

    # process noise
    w = get_noise(sigma_wp, sigma_wv)

    # state update
    x = A.dot(x) + B.dot(u) + w
    if np.linalg.norm(x[2:4]) > v_max:
        x[2:4] = v_max*x[2:4]/np.linalg.norm(x[2:4])
    x_sim.append(x)

    # log possible states at k = 0
    X0_arr.append(X0)

# simulate with moving horizon
for k in range(N, n_sim):

    # measurement
    v = get_noise(sigma_vp, sigma_vv)
    y = x + v
    y_sim.append(y)

    # MHE
    u_vec = u_sim[-N:]
    y_vec = y_sim[-N:]

    try:
        xhat, sol, X0 = mhe(X_Nm1, W, V, A, B, u_vec, y_vec, Qinv, Rinv, S, settings=settings)
    except ValueError as e:
        print(e)
        break

    xhat_sim.append(xhat)
    sol_times.append(sol.run_time)
    startup_times.append(sol.startup_time)
    iters.append(sol.iter)
    print(f'k = {k}, sol_time = {sol.run_time}, iter = {sol.iter}, |x-xhat|_2 = {np.linalg.norm(x-xhat)}')

    # control input
    if k < u_data.shape[0]:
        u = u_data[k,:]
    else:
        u = u_data[-1,:]
    u_sim.append(u)

    # process noise
    w = get_noise(sigma_wp, sigma_wv)

    # state update
    x = A.dot(x) + B.dot(u) + w
    if np.linalg.norm(x[2:4]) > v_max:
        x[2:4] = v_max*x[2:4]/np.linalg.norm(x[2:4])
    x_sim.append(x)

    # propagate possible states at k = N-1
    X_Nm1 = state_update(X_Nm1, A, B, W, V, u_vec[0], y_vec[0], S)

    # log possible states at k = 0
    X0_arr.append(X0)

# position mean squared error
e_pos_y = np.array([np.linalg.norm(y_sim[k-1][0:2] - x_sim[k][0:2]) for k in range(1, n_sim)])
rms_pos_y = np.sqrt(np.mean(e_pos_y**2))
e_pos_mhe = np.array([np.linalg.norm(xhat_sim[k-1][0:2] - x_sim[k][0:2]) for k in range(1, n_sim)])
rms_pos_mhe = np.sqrt(np.mean(e_pos_mhe**2))

e_vel_y = np.array([np.linalg.norm(y_sim[k-1][2:4] - x_sim[k][2:4]) for k in range(1, n_sim)])
rms_vel_y = np.sqrt(np.mean(e_vel_y**2))
e_vel_mhe = np.array([np.linalg.norm(xhat_sim[k-1][2:4] - x_sim[k][2:4]) for k in range(1, n_sim)])
rms_vel_mhe = np.sqrt(np.mean(e_vel_mhe**2))

print(f'Position RMS error from measurements: {rms_pos_y}')
print(f'Position RMS error from MHE: {rms_pos_mhe}')
print(f'Velocity RMS error from measurements: {rms_vel_y}')
print(f'Velocity RMS error from MHE: {rms_vel_mhe}')

# print solution and startup times
print(f'Avg. solution time [ms]: {1000.*np.mean(sol_times)}')
print(f'Avg. startup time [ms]: {1000.*np.mean(startup_times)}')
print(f'Avg. iterations: {np.mean(iters)}')
print(f'Max. solution time [ms]: {1000.*np.max(sol_times)}')
print(f'Max. startup time [ms]: {1000.*np.max(startup_times)}')
print(f'Max. iterations: {np.max(iters)}')
print(f'Min. solution time [ms]: {1000.*np.min(sol_times)}')
print(f'Min. startup time [ms]: {1000.*np.min(startup_times)}')
print(f'Min. iterations: {np.min(iters)}')

# plot settings
textwidth_pt = 10
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

# plot trajectory and position uncertainty sets
inches_per_pt = 1 / 72.27
figsize = (505.89 * inches_per_pt, 0.45*505.89 * inches_per_pt)  # Convert pt to inches

with plt.rc_context(rc_context):
    fig = plt.figure(constrained_layout=True, figsize=figsize)

    # position
    ax1 = fig.add_subplot(121)
    zono.plot(zono.project_onto_dims(X0_arr[0], [0,1]), ax1, color='y', alpha=0.1)
    for k in range(1, n_sim):
        zono.plot(zono.project_onto_dims(X0_arr[k], [0,1]), ax1, color='b', alpha=0.1)
        ax1.plot(x_sim[k][0], x_sim[k][1], 'xk')
        ax1.plot(y_sim[k-1][0], y_sim[k-1][1], '.g')
        ax1.plot(xhat_sim[k-1][0], xhat_sim[k-1][1], '.r')

    ax1.set_xlabel(r'$x$ [m]')
    ax1.set_ylabel(r'$y$ [m]')
    ax1.legend([r'Initial Set', r'Feasible Set', r'True', r'Measurement', r'MHE'])   
    ax1.axis('equal')
    ax1.grid(alpha=0.2)
    ax1.set_title(r'Position', fontsize=textwidth_pt)

    # velocity
    ax2 = fig.add_subplot(122)
    zono.plot(zono.project_onto_dims(X0_arr[0], [2,3]), ax2, color='y', alpha=0.1)
    for k in range(1, n_sim):
        zono.plot(zono.project_onto_dims(X0_arr[k], [2,3]), ax2, color='b', alpha=0.1)
        ax2.plot(x_sim[k][2], x_sim[k][3], 'xk')
        ax2.plot(y_sim[k-1][2], y_sim[k-1][3], '.g')
        ax2.plot(xhat_sim[k-1][2], xhat_sim[k-1][3], '.r')

    ax2.set_xlabel(r'$v_x$ [m/s]')
    ax2.set_ylabel(r'$v_y$ [m/s]') 
    ax2.axis('equal')
    ax2.grid(alpha=0.2)
    ax2.set_title(r'Velocity', fontsize=textwidth_pt)

    if is_latex_installed():
        plt.savefig('mhe_traj.pgf')
    
    plt.show()