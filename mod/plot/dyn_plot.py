"""
for dynamic plotting (third body orbit)
"""

import matplotlib.pyplot as plt

import mod.utils.static_utils as s_utils
import mod.plot.static_plot as s_plot
import mod.utils.dyn_utils as d_utils
import mod.consts as consts


def plot_dynamic_anim():
    """
    solve for the initial conditions given some jacobi integral
    """

    xy = d_utils.RK4([consts.RX0,consts.RY0,consts.VX0,consts.VY0])

    # plot the positions of the third body
    fig, ax = plt.subplots()

    s_plot.plot_lagrange_points()
    s_plot.plot_ZVC(s_utils.get_J("L1"))

    plt.plot(xy[0],xy[1])
    plt.show()