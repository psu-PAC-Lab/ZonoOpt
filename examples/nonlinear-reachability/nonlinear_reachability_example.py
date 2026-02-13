import zonoopt as zono
import numpy as np
import matplotlib.pyplot as plt

def f_unicycle(x, u, dt):
    """
    Discrete time unicycle dynamics
    x = [x, y, theta]
    u = [v, omega]
    """
    xp_k = x[0]
    yp_k = x[1]
    th_k = x[2]

    v_k = u[0]
    om_k = u[1]
    
    xp_kp1 = xp_k + v_k * np.cos(th_k) * dt
    yp_kp1 = yp_k + v_k * np.sin(th_k) * dt
    th_kp1 = th_k + om_k * dt
    return np.array([xp_kp1, yp_kp1, th_kp1])

def df_dxu_unicycle(x, u, dt):
    xp_k = x[0]
    yp_k = x[1]
    th_k = x[2]

    v_k = u[0]
    om_k = u[1]

    df_dxu = np.array([[1., 0., -v_k * np.sin(th_k) * dt, np.cos(th_k) * dt, 0.],
                      [0., 1.,  v_k * np.cos(th_k) * dt, np.sin(th_k) * dt, 0.],
                      [0., 0., 1., 0., dt]])
    return df_dxu

def df_dxu_unicycle_interval(x_int, u_int, dt):
    xp_k_int = x_int[0]
    yp_k_int = x_int[1]
    th_k_int = x_int[2]

    v_k_int = u_int[0]
    om_k_int = u_int[1]

    df_dxu_int = [[zono.Interval(1., 1.), zono.Interval(0., 0.), v_k_int * th_k_int.sin() * (-dt), th_k_int.cos() * dt, zono.Interval(0., 0.)],
                    [zono.Interval(0., 0.), zono.Interval(1., 1.), v_k_int * th_k_int.cos() * dt, th_k_int.sin() * dt, zono.Interval(0., 0.)],
                    [zono.Interval(0., 0.), zono.Interval(0., 0.), zono.Interval(1., 1.), zono.Interval(0., 0.), zono.Interval(dt, dt)]]
    df_dxu_lb = np.zeros((3, 5))
    df_dxu_ub = np.zeros((3, 5))
    for i in range(3):
        for j in range(5):
            df_dxu_lb[i,j] = df_dxu_int[i][j].lb
            df_dxu_ub[i,j] = df_dxu_int[i][j].ub
    return zono.IntervalMatrix(df_dxu_lb, df_dxu_ub)

def rand_x0(r_x0, r_th0):
    phi0 = np.random.uniform(0., 2*np.pi)
    r0 = np.random.uniform(0., r_x0)
    th0 = np.random.uniform(-r_th0, r_th0)
    xp0 = r0 * np.cos(phi0)
    yp0 = r0 * np.sin(phi0)
    return np.array([xp0, yp0, th0])

def rand_u(v_min, v_max, om_min, om_max):
    v = np.random.uniform(v_min, v_max)
    om = np.random.uniform(om_min, om_max)
    return np.array([v, om])

def monte_carlo_sim(n_samples, dt, N, r_x0, r_th0, v_min, v_max, om_min, om_max):
    x_monte_carlo = []
    x0_arr = []
    for _ in range(n_samples):
        x0_arr.append(rand_x0(r_x0, r_th0))
    x_monte_carlo.append(x0_arr)
    for k in range(N):
        xk_arr = []
        for xk in x_monte_carlo[-1]:
            u = rand_u(v_min, v_max, om_min, om_max)
            xkp1 = f_unicycle(xk, u, dt)
            xk_arr.append(xkp1)
        x_monte_carlo.append(xk_arr)
    return x_monte_carlo

# problem parameters
dt = 0.1
r_x0 = 0.01
r_th0 = 0.01
v_min = 1.
v_max = 2.
om_min = -0.1
om_max = 0.1
U = zono.interval_2_zono(zono.Box([v_min, om_min], [v_max, om_max]))
X0 = zono.cartesian_product(zono.make_regular_zono_2D(r_x0, 12, outer_approx=True), zono.interval_2_zono(zono.Box([-r_th0], [r_th0])))
N = 20

# compute reachable set
X_arr = [X0]
for k in range(N):
    Xk = X_arr[-1]
    f_xk = f_unicycle(Xk.get_center(), U.get_center(), dt)
    df_dxu_k = df_dxu_unicycle_interval(Xk.bounding_box(), U.bounding_box(), dt)
    
    XU = zono.cartesian_product(Xk, U)
    dXU = zono.minkowski_sum(XU, zono.Point(-XU.get_center()))
    Xkp1 = zono.affine_inclusion(dXU, df_dxu_k, f_xk)
    X_arr.append(Xkp1)

# plot reachable set and monte carlo
fig = plt.figure(constrained_layout=True)
ax = fig.add_subplot(111)
h = [None, None]
for X in X_arr:
    h[0] = zono.plot(zono.project_onto_dims(X, [0, 1]), ax=ax, color='b', alpha=0.1)[0]
for xk_arr in monte_carlo_sim(1000, dt, N, r_x0, r_th0, v_min, v_max, om_min, om_max):
    for x in xk_arr:
        h[1] = ax.plot(x[0], x[1], '.k', markersize=2)[0]
ax.axis('equal')
ax.grid(alpha=0.2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend(h, ['Bounding set', 'Monte Carlo samples'])

plt.show()
    