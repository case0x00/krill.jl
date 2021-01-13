"""
utils for the main krill script, supporting animating the
third body

"""

import math
from tqdm import tqdm
from numpy import linalg
import numpy as np

import mod.consts as consts

FUNC = "twobody"

if FUNC == "twobody":
    """
    two body simulation
    """

    TWOBODY_MU = 5.972e24 * 6.67408e-11 
    def f(rx, ry, vx, vy):
        """
        twobody equations of motion

        args:
            - state
        returns:
            - state
        """
        ax = -(TWOBODY_MU * rx) / ((rx ** 2 + ry ** 2) ** (3 / 2)) 
        ay = -(TWOBODY_MU * ry) / ((rx ** 2 + ry ** 2) ** (3 / 2)) 
        return [vx, vy, ax, ay]

elif FUNC == "CR3BP":
    """
    circular restricted three body problem simulation
    """


    def f(rx, ry, vx, vy):
        """
        the CR3BP equations of motion

        args:
            - state
        returns:
            - state
        """

        rb1 = [-consts.MU,0]
        rb2 = [1-consts.MU,0]

        r = [rx, ry]

        r1 = [a-b for a,b in zip(r, rb1)]
        r2 = [a-b for a,b in zip(r, rb2)]

        r13 = linalg.norm(r1,2)**3
        r23 = linalg.norm(r2,2)**3

        c1 = (1-consts.MU)/r13
        c2 = consts.MU/r23

        ax = 2 * vy + rx - c1 * (rx + consts.MU) - c2 * (rx - (1-consts.MU))
        ay = -2 * vx + ry - c1 * ry - c2 * ry

        return [vx, vy, ax, ay]

else:
    print("FUNC NOT IMPLEMENTED")



def RK4(rx, ry, vx, vy, h):
    """
    applies the RK4 to approximate solutions to the
    equations of motion
    """

    k1 = f(rx,ry,vx,vy)
    k2 = f(rx + (h * k1[0])/2, ry + (h * k1[1])/2, vx + (h * k1[2])/2, vy + (h * k1[3])/2)
    k3 = f(rx + (h * k2[0])/2, ry + (h * k2[1])/2, vx + (h * k2[2])/2, vy + (h * k2[2])/2)
    k4 = f(rx + h * k3[0], ry + h * k3[1], vx + h * k3[2], vy + h * k3[3])

    rx_step = rx + (h/6) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])    
    ry_step = ry + (h/6) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    vx_step = vx + (h/6) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    vy_step = vy + (h/6) * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])

    return rx_step, ry_step, vx_step, vy_step


def sim(init_conds, t0, h, tf):
    """
    simulate the orbit
    args:
        - init_conds -- [rx0,ry0,vx0,vy0]
        - t0 -- initial time
        - h -- time step
        - tf -- final time
    returns:
        - state
    """
    rx, ry, vx, vy = init_conds
    time_list = np.linspace(t0, tf, int(tf/h) + 1)

    rx_out, ry_out = np.zeros(len(time_list)), np.zeros(len(time_list))
    vx_out, vy_out = np.zeros(len(time_list)), np.zeros(len(time_list))

    for i in range(0, len(time_list)):
        rx, ry, vx, vy = RK4(rx, ry, vx, vy, h)
        rx_out[i], ry_out[i], vx_out[i], vy_out[i] = rx, ry, vx, vy
    return rx_out, ry_out, vx_out, vy_out
