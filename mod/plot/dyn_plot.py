"""
for dynamic plotting (third body orbit)
"""

import matplotlib.pyplot as plt
import math
import numpy as np

import mod.utils.static_utils as s_utils
import mod.plot.static_plot as s_plot
import mod.utils.dyn_utils as d_utils
import mod.consts as consts


def plot_dynamic_anim():
    """
    solve for the initial conditions given some jacobi integral
    """

    t0 = 0
    h = 10
#    tf = 24*60*60
    tf = 2

    x0 = [consts.RX0,consts.RY0,consts.VX0,consts.VY0]
    print(f"INITIAL CONDS: {x0}")
    rx, ry, _, _ = d_utils.sim(x0, t0, h, tf)

    # plot the positions of the third body
    fig, ax = plt.subplots()

    #s_plot.plot_lagrange_points()
    #s_plot.plot_ZVC(s_utils.get_J("L1"))

    plt.plot(rx,ry)
    plt.show()