"""
utils for the main krill script, supporting animating the
third body
"""

import numpy as np
from mod.utils.static import get_J_state
from mod.utils.solvers import RK4

import mod.consts as consts

FUNC = "CR3BP"

if FUNC == "twobody":
    """
    two body simulation
    """

    TWOBODY_MU = 5.972e24 * 6.67408e-11 
    def f(rx,ry,vx,vy):
        """
        twobody equations of motion

        args:
            - state
        returns:
            - state
        """
        ax = -(TWOBODY_MU * rx) / ((rx ** 2 + ry ** 2) ** (3 / 2)) 
        ay = -(TWOBODY_MU * ry) / ((rx ** 2 + ry ** 2) ** (3 / 2)) 
        return [vx,vy,ax,ay]

elif FUNC == "CR3BP":
    """
    circular restricted three body problem simulation
    """

    def f(rx,ry,vx,vy):
        """
        the CR3BP equations of motion

        args:
            - state
        returns:
            - state
        """

        r1 = ((rx+consts.MU)**2 + ry**2)**(1/2)
        r2 = ((rx-(1-consts.MU))**2 + ry**2)**(1/2)

        m1 = 1-consts.MU
        m2 = consts.MU

        ax = rx + 2*vy + m1*(-consts.MU-rx)/(r1**3) + m2*(1-consts.MU-rx)/(r2**3)
        ay = ry - 2*vx - m1*ry/(r1**3) - m2*ry/(r2**3)
        
        # x^2 + y^2 + 2(1-mu)/r1 + 2mu/r2 - (xdot^2 + ydot^2) = C


        return [vx,vy,ax,ay]


    def neg_f(rx,ry,vx,vy):
        """
        reverse equations of motion
        """
        return [-1*i for i in f(rx,ry,vx,vy)]


else:
    print("FUNC NOT IMPLEMENTED")




def sim(init_conds, t0, h, tf):
    """
    simulate the orbit

    args:
        - init_conds: [rx0,ry0,vx0,vy0]
        - t0: initial time
        - h: time step
        - tf: final time
    returns:
        - state
    """
    rx, ry, vx, vy = init_conds
    time_list = np.linspace(t0, tf, int(tf/h) + 1)

    rx_out, ry_out = np.zeros(len(time_list)), np.zeros(len(time_list))
    vx_out, vy_out = np.zeros(len(time_list)), np.zeros(len(time_list))

    for i in range(0, len(time_list)):
        rx, ry, vx, vy = RK4(f, rx, ry, vx, vy, h, i)
        print(f"J: {get_J_state([rx,ry,vx,vy])}")
        rx_out[i], ry_out[i], vx_out[i], vy_out[i] = rx, ry, vx, vy
    return rx_out, ry_out, vx_out, vy_out
