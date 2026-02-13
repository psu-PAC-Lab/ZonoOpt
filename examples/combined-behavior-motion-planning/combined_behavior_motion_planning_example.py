import zonoopt as zono
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.signal import cont2discrete as c2d
from enum import Enum
import copy
import control
from draw_car import draw_car

# globals:
MODE = 'plot_scenario' # options: 'warmstart_test', 'plot_scenario'
CAR_LENGTH = 0.2
TEXTWIDTH_PT = 10
RC_CONTEXT = {
    "text.usetex": True,
    "font.size": TEXTWIDTH_PT,
    "font.family": "serif",  # Choose a serif font like 'Times New Roman' or 'Computer Modern'
    "pgf.texsystem": "pdflatex",
    "pgf.rcfonts": False,
}
INCHES_PER_POINT = 1 / 72.27


class Lane(Enum):
    LEFT = 0
    RIGHT = 1

class Car():
    def __init__(self, length, position, lane, vel):
        """
        Initialize Car object

        Args:
            length (float): length of the car for purposes of planning (includes buffer)
            position (float): position of the car along the road
            lane (Lane): lane of the car (Lane.LEFT or Lane.RIGHT)
            vel (float): velocity of car
        """
        self.length = length
        self.position = position
        self.lane = lane
        self.vel = vel

    def __copy__(self):
        return Car(self.length, self.position, self.lane, self.vel)

    def __deepcopy__(self, memo):
        if self in memo:
            return memo[self]
        new_car = Car(self.length, self.position, self.lane, self.vel)
        memo[self] = new_car

        return new_car
    
    def __repr__(self):
        if self.lane == Lane.LEFT:
            lane_str = 'left'
        else:
            lane_str = 'right'
        return f'Car: length = {self.length}, position = {self.position}, lane = {lane_str}, vel = {self.vel}'


class TrafficScenario():
    def __init__(self, cars, n_time_steps, dt, min_zone_length, road_length, lane_width):
        
        # fields
        self.n_time_steps = n_time_steps
        self.dt = dt
        self.y_ll = -lane_width/2. 
        self.y_rl = lane_width/2.
        self.road_start = 0.
        self.road_length = road_length
        self.lane_width = lane_width
        self.min_zone_length = min_zone_length

        # init
        cars_k = copy.deepcopy(cars)
        self.cars_arr = []
        self.zones_arr = []
        self.hz_arr = []
        self.cars_hz_arr = []

        for _ in range(n_time_steps+1):
        
            # get car positions and append
            self.cars_arr.append(copy.deepcopy(cars_k))

            # get zones
            ll_cars = [car for car in cars_k if car.lane == Lane.LEFT]
            rl_cars = [car for car in cars_k if car.lane == Lane.RIGHT]
            self.zones_arr.append([])
            self._get_zones(self.zones_arr[-1], ll_cars, Lane.LEFT)
            self._get_zones(self.zones_arr[-1], rl_cars, Lane.RIGHT)

            # create zones as HZs
            self.hz_arr.append(self._zones_2_hz(self.zones_arr[-1]))

            # create cars as HZs
            car_zonos = []
            for car in cars_k:
                car_zonos.append(self._car_to_z(car, lane_width))
            if len(car_zonos) == 0:
                self.cars_hz_arr.append(zono.EmptySet(2))
            else:
                self.cars_hz_arr.append(zono.zono_union_2_hybzono(car_zonos))

            # update car positions
            for car in cars_k:
                car.position += car.vel * dt

    def get_hzs(self, k):
        assert k >= 0 and k <= self.n_time_steps
        return self.hz_arr[k]
    
    def get_cars(self, k):
        assert k >= 0 and k <= self.n_time_steps
        return self.cars_arr[k]
    
    def get_free_space(self, k):
        return self.get_hzs(k)
    
    def get_all_space(self):
        return zono.interval_2_zono(zono.Box([self.road_start, -self.lane_width], [self.road_length, self.lane_width]))

    @staticmethod
    def make_random(n_cars, car_length, car_speed_mean, car_speed_std, n_time_steps, dt, min_zone_length, road_length, lane_width, seed=0):
        cars = TrafficScenario._generate_cars(n_cars=n_cars, 
                                              car_length=car_length, 
                                              car_speed_mean=car_speed_mean, 
                                              car_speed_std=car_speed_std, 
                                              road_length=road_length, 
                                              n_time_steps=n_time_steps, 
                                              dt=dt, 
                                              seed=seed)
        return TrafficScenario(cars=cars, 
                               n_time_steps=n_time_steps, 
                               dt=dt, 
                               min_zone_length=min_zone_length, 
                               road_length=road_length, 
                               lane_width=lane_width)

    def plot(self, k, ax):
        assert k>=0 and k<=self.n_time_steps
        zono.plot(self.hz_arr[k], ax=ax, color='b', alpha=0.1, edgecolor='k')
        
    @staticmethod
    def _generate_cars(n_cars, car_length, car_speed_mean, car_speed_std, road_length, n_time_steps, dt, seed):
        np.random.seed(seed)

        cars = []
        while len(cars) < n_cars:
            new_car = Car(length=car_length, 
                          position=np.random.uniform(car_length/2., road_length - car_length/2.), 
                          lane=Lane.LEFT if np.random.rand() < 0.5 else Lane.RIGHT, 
                          vel=np.random.normal(car_speed_mean, car_speed_std))
            if TrafficScenario._check_cars_valid(cars, new_car, n_time_steps, dt):
                cars.append(new_car)
                # print(new_car)
        return cars

    def _get_zones(self, zones, cars, lane):

        # sort cars
        cars.sort(key=lambda car: car.position)

        # find zones
        for i in range(len(cars)):
            end = cars[i].position - cars[i].length/2.
            end = min(max(end, self.road_start), self.road_length)

            start = self.road_start
            if i > 0:
                start = max(start, cars[i-1].position + cars[i-1].length/2.)

            if end - start > self.min_zone_length:
                zones.append((start, end, lane))

        # after last car
        if len(cars) == 0:
            start = self.road_start
        else:
            start = max(cars[-1].position + cars[-1].length/2., self.road_start)
        end = self.road_length
        if end - start > self.min_zone_length:
            zones.append((start, end, lane))      

    def _zones_2_hz(self, zones):
        zonos = []
        for zone in zones:
            if zone[2] == Lane.LEFT:
                lane_min = -self.lane_width
                lane_max = 0.
            else:
                lane_min = 0.
                lane_max = self.lane_width
            zonos.append(zono.interval_2_zono(zono.Box([zone[0], lane_min], [zone[1], lane_max])))
        if len(zonos) == 0:
            return zono.EmptySet(2)
        else:
            return zono.zono_union_2_hybzono(zonos)
        
    def _car_to_z(self, car, lane_width):
        if car.lane == Lane.LEFT:
            lane_min = -lane_width
            lane_max = 0.
        else:
            lane_min = 0.
            lane_max = lane_width
        return zono.interval_2_zono(zono.Box([car.position - car.length/2., lane_min], [car.position + car.length/2., lane_max]))
    
    @staticmethod
    def _check_cars_valid(cars, new_car, n_time_steps, dt):
        cars_cp = copy.deepcopy(cars)
        new_car_cp = copy.copy(new_car)

        for k in range(n_time_steps+1):

            # check conflicts
            for car in cars_cp:
                if car.lane == new_car_cp.lane and abs(car.position - new_car_cp.position) < (car.length/2. + new_car_cp.length/2.):
                    return False
            
            # update car positions
            for car in cars_cp:
                car.position += car.vel * dt
            new_car_cp.position += new_car_cp.vel * dt
                
        return True


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

def motion_planning_problem(settings=zono.OptSettings(), rng_seed=None, plot_constraints=False):
    
    # traffic scenario

    # speed of traffic
    vref = 0.5

    # discrete time step
    dt = 1.0

    # number of time steps
    N = 15

    # velocity and heading bounds
    vmax = 0.7

    # max lane deviation angle
    mu_max = 30.*(np.pi/180.)

    # control bounds
    a_bounds = (-1., 1.)
    dddot_pert_bounds = (-0.01, 0.01)

    # scenario
    car_length_buf = 0.2 + 2. * vmax * dt
    road_length = vmax * N * dt

    car_rl = Car(length=car_length_buf, position=3.5, lane=Lane.RIGHT, vel=0.05)
    car_ll = Car(length=car_length_buf, position=2.5, lane=Lane.LEFT, vel=0.3)
    cars = (car_rl, car_ll)
    
    if rng_seed is None:
        scenario = TrafficScenario(cars=cars, 
                                n_time_steps=N, 
                                dt=dt,
                                min_zone_length=0.1*car_length_buf, 
                                road_length=road_length, 
                                lane_width=0.51)
    else:
        scenario = TrafficScenario.make_random(2, car_length_buf, 0.2, 0.1, N, dt, 0.1*car_length_buf, road_length, 0.51, seed=rng_seed)

    # ego car
    ego_car = Car(length=0.2, position=0.0, lane=Lane.RIGHT, vel=vref)
    
    # make sure that ego car is at least not inside another car within first 3 sec
    k3 = int(np.ceil(3./dt))
    for k in range(k3):
        cars_init = scenario.get_cars(k)
        for car in cars_init:
            if car.lane == ego_car.lane and abs(car.position - ego_car.position) < (car.length/2. + ego_car.length/2.):
                raise RuntimeError(f'Invalid planning problem: ego car is in collision with another car at k={k}')

    # dynamics model

    # lateral dynamics in Frenet frame using differential flatness
    # xlat = [d, ddot], ulat = [dddot]
    A_lat = np.array([[1., dt],
                      [0., 1.]])
    B_lat = np.array([[0.5*dt**2],
                      [dt]])
    
    # lane-following control law: dddot = -[kd, kddot]*[d-dr; ddot]
    Q_lat = np.diag([1., 0.])
    R_lat = np.array([[10.]])
    K_lat, _, _ = control.dlqr(A_lat, B_lat, Q_lat, R_lat)
    kd = K_lat[0,0]
    kddot = K_lat[0,1]
    # print(f'kd = {kd}, kddot = {kddot}')

    # simulate to get max ddot during a turn
    x = np.array([scenario.y_rl, 0.]) # start in right lane
    xr = np.array([scenario.y_ll, 0.]) # target left lane
    ddot_max = 0. # init
    while np.abs(x[0] - xr[0]) > 0.01 * np.abs(scenario.y_rl - scenario.y_ll):
        u = -K_lat.dot(x - xr)
        x = A_lat.dot(x) + B_lat.dot(u)
        ddot_max = max(ddot_max, np.abs(x[1]))

    # print(f'Max ddot during lane change: {ddot_max} m/s2')

    # min sdot for turning
    sdot_min_turn = ddot_max / np.tan(mu_max)
    # print(f'Min sdot for turning: {sdot_min_turn} m/s')

    # full dynamics model
    # x = [s, d, sdot, ddot]
    # u = [sddot, dddot]
    # x_k+1 = A*x_k + B*u_k
    A = np.array([[1., 0., dt, 0.],
                    [0., 1., 0., dt],
                    [0., 0., 1., 0.],
                    [0., 0., 0., 1.]])
    B = np.array([[0.5*dt**2, 0.],
                    [0., 0.5*dt**2],
                    [dt, 0.],
                    [0., dt]])

    # closed-loop dynamics
    K = np.array([[0., 0., 0., 0.],
                  [0., kd, 0., kddot]])
    A_rl = A - B @ K
    B_rl = B
    E_rl = (B @ K @ np.array([[0.], [scenario.y_rl], [0.], [0.]])).flatten()

    A_ll = A - B @ K
    B_ll = B
    E_ll = (B @ K @ np.array([[0.], [scenario.y_ll], [0.], [0.]])).flatten()

    # get state and input constraints
    eps_ddot = 0.3 * ddot_max
    ddot_p = ((ddot_max - eps_ddot) / sdot_min_turn) * vmax + eps_ddot

    S = zono.cartesian_product(scenario.get_all_space(),
        zono.interval_2_zono(zono.Box([0., -vmax], [vmax, vmax])))
    
    S_ll = zono.cartesian_product(scenario.get_all_space(),
        zono.vrep_2_conzono(np.array([[0., -eps_ddot],
                                      [0., 0.],
                                      [vmax, -ddot_p],
                                      [vmax, 0.]])))
    S_rl = zono.cartesian_product(scenario.get_all_space(),
        zono.vrep_2_conzono(np.array([[0., 0.],
                                      [0., eps_ddot],
                                      [vmax, 0.],
                                      [vmax, ddot_p]])))
    
    if plot_constraints:
        fig_width_pt = 245.71
        figsize = (fig_width_pt * INCHES_PER_POINT, 0.6 * fig_width_pt * INCHES_PER_POINT)

        with plt.rc_context(RC_CONTEXT):
            fig = plt.figure(constrained_layout=True, figsize=figsize)  
            ax = fig.add_subplot(1,1,1)

            h = []
            h.append(zono.plot(zono.project_onto_dims(S_rl, [2,3]), ax=ax, color='r', alpha=0.2, edgecolor='k')[0])
            h.append(zono.plot(zono.project_onto_dims(S_ll, [2,3]), ax=ax, color='b', alpha=0.2, edgecolor='k')[0])
            ax.set_xlabel(r'$\dot{s}$~[m/s]')
            ax.set_ylabel(r'$\dot{d}$~[m/s]')
            ax.grid(alpha=0.2)
            ax.legend(h, (r'$S^{\mathit{rl}}$', r'$S^{\mathit{ll}}$'), loc='upper left', fontsize=TEXTWIDTH_PT)

            plt.savefig('behavior_motion_planning_state_constraints.pgf')
            plt.show()

    U = zono.interval_2_zono(zono.Box([a_bounds[0], dddot_pert_bounds[0]], [a_bounds[1], dddot_pert_bounds[1]]))

    # get graph of function for all time steps
    Psi_arr = []
    nx = 4
    nu = 2
    for k in range(N):
        Psi_rl = make_graph_of_function(A_rl, B_rl, E_rl, zono.cartesian_product(S_rl, U))
        Psi_ll = make_graph_of_function(A_ll, B_ll, E_ll, zono.cartesian_product(S_ll, U))
        if Psi_rl.is_zono() and Psi_ll.is_zono():
            Psi = zono.zono_union_2_hybzono((Psi_rl, Psi_ll))
        else:
            Psi = zono.union_of_many((Psi_rl, Psi_ll), preserve_sharpness=False)
        Psi = zono.intersection_over_dims(Psi, scenario.get_free_space(k+1), [i for i in range(nx+nu, nx+nu+2)])
        
        Psi_arr.append(Psi)

    # build and solve planning problem
    if ego_car.lane == Lane.LEFT:    
        x0 = np.array([ego_car.position, scenario.y_ll, 0., 0.])
    else:
        x0 = np.array([ego_car.position, scenario.y_rl, 0., 0.])
    
    xr_arr = []
    s = ego_car.position
    for _ in range(N):
        s += vref * dt
        xr_arr.append(np.array([s, scenario.y_rl, vref, 0.]))

    Q = np.diag([0., 0., 0.5, 0.5])
    R = np.diag([10., 10.])
    QN = np.diag([10., 10., 10., 10.])

    sol = zono.OptSolution()
    x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u = zono_planning_prob(x0, xr_arr, Psi_arr, S, U, Q, R, QN, N, settings=settings, sol=sol)

    return x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u, sol, scenario

def warmstart_plan(Z_prob, z_ws, P, q, c, idx_x, idx_u, settings=zono.OptSettings(), rng_seed=0, pert_mag=0.1):
    """
    Warmstart planning problem

    Args:
        Z_prob (zono.HybZono): planning problem set
        P (scipy.sparse.csc_matrix): state cost matrix
        q (np.ndarray): state cost vector
        c (float): cost offset
        settings (zono.OptSettings, optional): optimization settings
        sol (zono.OptSolution, optional): optimization solution reference

    Returns:
        ...
    """

    # perturb z_ws to test robustness
    np.random.seed(rng_seed)
    z_ws += pert_mag * np.random.randn(z_ws.shape[0])

    # convex relaxation
    Z_prob_cvx = Z_prob.convex_relaxation()
    sol_cvx = zono.OptSolution()
    Z_prob_cvx.project_point(z_ws, solution=sol_cvx)

    # warmstart and resolve
    sol = zono.OptSolution()
    warm_start_params = zono.WarmStartParams()
    warm_start_params.z = sol_cvx.z
    warm_start_params.u = sol_cvx.u
    z_opt_ws = Z_prob.optimize_over(P, q, c=c, settings=settings, solution=sol, warm_start_params=warm_start_params)

    # state, input trajectory
    x_traj = []
    u_traj = []
    for idx in idx_x:
        x_traj.append(z_opt_ws[idx])
    for idx in idx_u:
        u_traj.append(z_opt_ws[idx])

    # include warm-start time in total runtime
    sol.run_time += sol_cvx.run_time

    return x_traj, u_traj, z_opt_ws, sol


def plot_scenario(scenario, x_traj, n_subplots, n_cols=1, height_ratio=0.25):

    figsize = (505.89 * INCHES_PER_POINT, height_ratio*505.89 * INCHES_PER_POINT)
    
    with plt.rc_context(RC_CONTEXT):

        fig = plt.figure(constrained_layout=True, figsize=figsize)
        
        n_rows = int(np.ceil(n_subplots / n_cols))
        gs = fig.add_gridspec(n_rows, n_cols)

        for i in range(n_rows):
            for j in range(n_cols):

                cnt = i * n_cols + j
                if cnt >= n_subplots:
                    continue

                ax = fig.add_subplot(gs[i,j])

                ax.invert_yaxis()
                k_plt = int(np.floor(cnt/(n_subplots-1) * (scenario.n_time_steps)))
                scenario.plot(k=k_plt, ax=ax)

                for k in range(k_plt+1):
                    # obstacle cars
                    for car in scenario.get_cars(k):

                        
                        if car.lane == Lane.LEFT:
                            d = scenario.y_ll
                        else:
                            d = scenario.y_rl
                        s = car.position
                        th = 0.
                        if k == k_plt:
                            alpha = 1.
                        else:
                            alpha = 0.25

                        draw_car(ax, s, d, th, 0., CAR_LENGTH, scenario.lane_width/4., scenario.lane_width/6., 'red', alpha=alpha)
 
                    # ego car
                    s = x_traj[k][0]
                    d = x_traj[k][1]
                    th = np.arctan2(x_traj[k][3], x_traj[k][2])
                    lanewidth = scenario.lane_width
                    if k == k_plt:
                        alpha = 1.
                    else:
                        alpha = 0.25

                    draw_car(ax, s, d, th, 0., CAR_LENGTH, lanewidth/4., lanewidth/6., 'green', alpha=alpha)

                ax.set_title(rf'$k = {k_plt}$', fontsize=TEXTWIDTH_PT)
                ax.axis('equal')
                ax.grid(alpha=0.2)

                if i < n_rows - 1:
                    ax.set_xticklabels([])
                if j != 0:
                    ax.set_yticklabels([])

        fig.supxlabel(r'$s$ [m]', fontsize=TEXTWIDTH_PT)
        fig.supylabel(r'$d$ [m]', fontsize=TEXTWIDTH_PT)

        plt.savefig('behavior_motion_planning_scenario.pgf')
        plt.show()

def plot_warmstart_stats(cr_iter_arr, cr_time_arr, ws_iter_arr, ws_time_arr):
    
    figwidth_pt = 245.71
    figheight_pt = 0.5*figwidth_pt
    
    with plt.rc_context(RC_CONTEXT):
        fig = plt.figure(constrained_layout=True, figsize=(figwidth_pt * INCHES_PER_POINT, figheight_pt * INCHES_PER_POINT))

        colors = ('#1f77b4', '#ff7f0e')

        # iterations
        ax = fig.add_subplot(1,2,1)
        bp = ax.boxplot(
            [cr_iter_arr, ws_iter_arr],
            tick_labels=[r'C.R.', r'W.S.'],
            patch_artist=True,
            medianprops={'color': 'black'},
            showfliers=True,
            sym='k.',
            whis=(0., 100.) # whiskers cover all data)
        )
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.5)

        ax.set_title(r'Iterations', fontsize=TEXTWIDTH_PT)
        ax.set_ylabel(r'[count]')
        ax.set_yscale('log')
        ax.grid(axis='y', which='major', alpha=0.2)

        # run time
        ax = fig.add_subplot(1,2,2)
        bp = ax.boxplot(
            [cr_time_arr, ws_time_arr],
            tick_labels=[r'C.R.', r'W.S.'],
            patch_artist=True,
            medianprops={'color': 'black'},
            showfliers=True,
            sym='k.',
            whis=(0., 100.) # whiskers cover all data)
        )
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.5)

        ax.set_title(r'Run Time', fontsize=TEXTWIDTH_PT)
        ax.set_ylabel(r'[sec]')
        ax.set_yscale('log')
        ax.grid(axis='y', which='major', alpha=0.2)
        labels = ax.get_xticklabels()
        for label in labels:
            label.set_horizontalalignment('center')

        plt.savefig('behavior_motion_planning_warmstart_stats.pgf')
        plt.show()


### MAIN ###

# settings
settings = zono.OptSettings()
settings.t_max = 1.
settings.k_max_admm_fp_ph1 = 5000
settings.k_max_admm_fp_ph2 = 1000000
settings.verbosity_interval = 1000
settings.single_threaded_admm_fp = True
settings.polish = False
settings.enable_rng_seed = True
settings.rng_seed = 0
settings.eps_prim = 0.01
settings.eps_dual = 0.1
settings.k_restart = 1000
settings.rho = 100.

if MODE == 'warmstart_test':

    settings.verbose = False

    n_trials = 100
    trial = 0
    seed = 0
    cr_iter_arr = []
    cr_time_arr = []
    cr_obj_arr = []
    ws_iter_arr = []
    ws_time_arr = []
    ws_obj_arr = []
    cr_cnt = 0
    ws_cnt = 0

    while trial < n_trials:

        # increment seed
        settings.rng_seed = seed
        seed += 1

        # build and solve problem
        try:
            x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u, sol_cr, scenario = motion_planning_problem(settings, rng_seed=seed)    
        except RuntimeError:
            continue

        print(f'Trial: {trial}')
        trial += 1
        
        if sol_cr.converged:
            cr_cnt += 1
            cr_iter_arr.append(sol_cr.iter)
            cr_time_arr.append(sol_cr.run_time)
            cr_obj_arr.append(sol_cr.J)
        else:
            continue

        # warmstart and retry
        x_traj, u_traj, z_opt_ws, sol_ws = warmstart_plan(Z_prob, z_opt, P, q, c, idx_x, idx_u, settings=settings, rng_seed=seed)
        if sol_ws.converged:
            ws_cnt += 1
            ws_iter_arr.append(sol_ws.iter)
            ws_time_arr.append(sol_ws.run_time)        
            ws_obj_arr.append(sol_ws.J) 

    # print stats
    print(f'Convex relaxation: iters: median={np.median(cr_iter_arr)} s, time: median={np.median(cr_time_arr)} s')
    print(f'Warmstart: iters: median={np.median(ws_iter_arr)} s, time: median={np.median(ws_time_arr)} s')
    print(f'Convex relaxation converged in {cr_cnt}/{n_trials} trials')
    print(f'Warmstart converged in {ws_cnt}/{cr_cnt} trials')
    print(f'Objectives: CR median={np.median(cr_obj_arr)}, WS median={np.median(ws_obj_arr)}')

    # plot stats
    plot_warmstart_stats(cr_iter_arr, cr_time_arr, ws_iter_arr, ws_time_arr)
    
elif MODE == 'plot_scenario':
    settings.verbose = True

    x_traj, u_traj, z_opt, Z_prob, P, q, c, idx_x, idx_u, sol_cr, scenario = motion_planning_problem(settings, plot_constraints=False)
    plot_scenario(scenario, x_traj, 6, n_cols=2, height_ratio=0.35)

else:
    raise RuntimeError(f'Invalid mode: {MODE}')