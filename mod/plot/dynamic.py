"""
for dynamic plotting (third body orbit)
"""

import matplotlib.pyplot as plt
import math
import numpy as np

import mod.utils.static as s_utils
import mod.plot.static as s_plot
import mod.utils.dynamic as d_utils
import mod.consts as consts


def plot_anim():
    """
    solve for the initial conditions given some jacobi integral
    """

    t0 = 0
    h = 10
    tf = 24*60*60
    J = 3.18


    x0 = [consts.RX0,consts.RY0,consts.VX0,consts.VY0]
    print(f"INITIAL CONDS: {x0}")
    rx, ry, _, _ = d_utils.sim(x0, t0, h, tf)

    # plot the positions of the third body
    fig, ax = plt.subplots()
    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.grid("on",linestyle=":")
    ax.set_aspect("equal")
    

    s_plot.plot_libration_points()
    #print(s_utils.get_J_state(x0))

    s_plot.plot_ZVC(3.18)

    plt.plot(rx,ry)
    plt.show()