import matplotlib.pyplot as plt
import numpy as np

def draw_car(ax, x, y, th, psi, length, width, chassis_width, color, alpha=1.0):
    """
    Draws depiction of a car on given axis.

    Args:
        ax (matplotlib.axes.Axes): axis to draw on
        x (float): x position of car center
        y (float): y position of car center
        th (float): heading angle of car body
        psi (float): steering angle of front wheels
        length (float): total car length
        width (float): total car width
        chassis_width (float): width of car chassis (without wheels)
        color (str or tuple): color of car chassis
        alpha (float): transparency of fills

    Returns:
        tuple: (chassis_fill, wheel_fills, axle_fills)
            chassis_fill (list of matplotlib.patches.Polygon): fill for chassis
            wheel_fills (list of list of matplotlib.patches.Polygon): fills for wheels
            axle_fills (list of list of matplotlib.patches.Polygon): fills for axles
    """

    # wheel
    col_wheel = 'black'
    wheel_wid = 0.5 * width
    wheel_len = (1/3) * length

    # axle
    col_axle = (194/255, 197/255, 204/255)
    axle_len = (width - wheel_wid - chassis_width) / 2.
    axle_wid = wheel_len / 4.

    # rotation matrix
    C_th = np.array([[np.cos(-th), np.sin(-th)],
                    [-np.sin(-th), np.cos(-th)]])
    C_psi = np.array([[np.cos(-psi), np.sin(-psi)],
                    [-np.sin(-psi), np.cos(-psi)]])
    
    # chassis

    # geometry
    x_chassis = np.array([-length/2, length/2, length/2, -length/2])
    y_chassis = np.array([-width/2, -width/2, width/2, width/2])

    # rotate
    xy_chassis = C_th @ np.vstack((x_chassis, y_chassis))
    x_chassis = xy_chassis[0,:] + x
    y_chassis = xy_chassis[1,:] + y

    # make fill
    chassis_fill = ax.fill(x_chassis, y_chassis, color=color, edgecolor=None, alpha=alpha)

    # wheels

    # geometry
    x_wheel = np.array([-wheel_len/2, wheel_len/2, wheel_len/2, -wheel_len/2])
    y_wheel = np.array([-wheel_wid/2, -wheel_wid/2, wheel_wid/2, wheel_wid/2])

    # shifts of wheels relative to center
    wheel_name = ('RL', 'FL', 'FR', 'RR')  # rear-left, forward-right, etc.
    x_shift = np.array([-length/2 + wheel_len/2,
                        length/2 - wheel_len/2,
                        length/2 - wheel_len/2,
                        -length/2 + wheel_len/2])
    y_shift = np.array([-width/2 - axle_len - wheel_wid/2,
                        -width/2 - axle_len - wheel_wid/2,
                        width/2 + axle_len + wheel_wid/2,
                        width/2 + axle_len + wheel_wid/2])
    
    # get fills for each wheel
    x_wheels = np.zeros((4,4))
    y_wheels = np.zeros((4,4))
    wheel_fills = []
    for i in range(4):
        # rotated wheel in body frame
        if wheel_name[i] in ('FL', 'FR'):
            xy_wheel_body_rot = C_psi @ np.vstack((x_wheel, y_wheel))
        else:
            xy_wheel_body_rot = np.vstack((x_wheel, y_wheel))

        # rotate body
        xy_wheel_rot = C_th @ (np.array([[x_shift[i]], [y_shift[i]]]) + xy_wheel_body_rot)

        # offset by position of body
        x_wheels[i,:] = xy_wheel_rot[0,:] + x
        y_wheels[i,:] = xy_wheel_rot[1,:] + y

        # make fill
        wheel_fills.append(ax.fill(x_wheels[i,:], y_wheels[i,:], color=col_wheel, edgecolor=None, alpha=alpha))

    # axles

    # geometry
    x_axle = np.array([-axle_wid/2, axle_wid/2, axle_wid/2, -axle_wid/2])
    y_axle = np.array([-axle_len/2, -axle_len/2, axle_len/2, axle_len/2])

    # shifts of axles relative to center
    x_shift = np.array([-length/2 + wheel_len/2,
                        length/2 - wheel_len/2,
                        length/2 - wheel_len/2,
                        -length/2 + wheel_len/2])
    y_shift = np.array([-width/2 - axle_len/2,
                        -width/2 - axle_len/2,
                        width/2 + axle_len/2,
                        width/2 + axle_len/2])
    
    # get fills for each axle
    x_axles = np.zeros((4,4))
    y_axles = np.zeros((4,4))
    axle_fills = []
    for i in range(4):
        # rotate body
        xy_axle_body = np.vstack((x_axle, y_axle))
        xy_axle_rot = C_th @ (np.array([[x_shift[i]], [y_shift[i]]]) + xy_axle_body)

        # offset by position of body
        x_axles[i,:] = xy_axle_rot[0,:] + x
        y_axles[i,:] = xy_axle_rot[1,:] + y

        # make fill
        axle_fills.append(ax.fill(x_axles[i,:], y_axles[i,:], color=col_axle, edgecolor=None, alpha=alpha))

    return chassis_fill, wheel_fills, axle_fills