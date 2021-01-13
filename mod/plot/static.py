"""
for static plots (ZVC)
"""


import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math

import mod.consts as consts
import mod.utils.static as s_utils


def plot_libration_points():
    """
    computes the position of the 5 lagrange points
    https://www.mat.univie.ac.at/~westra/lagrangepoints.pdf
    """

    # L1 -- between the two masses
    # r1 = consts.MU, r2 = 1-consts.MU
    # x = a(1 - cbrt(M2/3M1)) -> x = (1-cbrt(consts.MU/3))
    L1 = [(1-math.pow(consts.MU/3,1/3)),0]
    #print(f"L1 at {L1}")
    plt.text(L1[0]-0.2,L1[1],"L1")
    plt.scatter(L1[0],L1[1],s=5,c="black")

    # L2 -- on the line of the masses, beyond the smaller one
    # x = r2(1 + cbrt(M2/3M1)) -> x = (1-consts.MU)(1+cbrt(consts.MU/3)) 
    L2 = [(1-consts.MU)*(1+math.pow(consts.MU/3,1/3)),0]
    #print(f"L2 at {L2}")
    plt.text(L2[0]+0.05,L2[1],"L2")
    plt.scatter(L2[0],L2[1],s=5,c="black")

    # L3 -- on the line of masses, beyond the larger one
    # x = r2(1 + 17/12 M2/M1) -> (1-consts.MU)(1 + 17/12 consts.MU)
    L3 = [-(1-consts.MU)*(1+17/12 * consts.MU),0]
    #print(f"L3 at {L3}")
    plt.text(L3[0]+0.05,L3[1],"L3")
    plt.scatter(L3[0],L3[1],s=5,c="black")

    # L4 -- corner of top equilateral triangle
    # eq 2.45 
    # https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf
    # x = 1/2 - consts.MU, y = sqrt3/2
    L4 = [1/2-consts.MU,math.sqrt(3)/2]
    #print(f"L4 at {L4}")
    plt.text(L4[0]+0.05,L4[1]+0.05,"L4")
    plt.scatter(L4[0],L4[1],s=5,c="black")

    # L5 -- corner of bottom equilateral triangle
    L5 = [1/2-consts.MU,-math.sqrt(3)/2]
    #print(f"L5 at {L5}")
    plt.text(L5[0]+0.05,L5[1]-0.13,"L5")
    plt.scatter(L5[0],L5[1],s=5,c="black")

def plot_ZVC(J):
    """
    plots the hill region (zero velocity curves).

    https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf
    eq 2.53

    args:
        J -- jacobi integral

    """

    x = np.linspace(consts.LOWERLIM, consts.UPPERLIM, 150)
    y = np.linspace(consts.LOWERLIM, consts.UPPERLIM, 150)

    X, Y = np.meshgrid(x,y)

    r1 = ((X+consts.MU)**2 + Y**2)**(1/2)
    r2 = ((X-(1-consts.MU))**2 + Y**2)**(1/2)

    Z = X**2 + Y**2 + (2*(1-consts.MU))/r1 + (2*consts.MU)/r2 - J

    plt.contour(X,Y,Z,[0],colors="black")
    
def plot():
    """
    builds 4 subplots with variations in jacobi integral
    to demonstrate the ZVC
    """

    fig = plt.figure(figsize=(10,10))

    ax = fig.add_subplot(221)

    # libration points
    plot_libration_points()

    # ZVC
    point = "L1"
    plot_ZVC(3.2)


    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    ax.set_aspect("equal")
    plt.title("$J > J_{L1}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(222)

    # libration points
    plot_libration_points()

    # ZVC
    point = "L2"
    plot_ZVC(3.180)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    ax.set_aspect("equal")
    plt.title("$J_{L2} < J < J_{L1}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(223)

    # libration points
    plot_libration_points()

    # ZVC
    point = "L3"
    plot_ZVC(3.10)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    ax.set_aspect("equal")
    plt.title("$J_{L3} < J < J_{L2}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(224)

    #libration points
    plot_libration_points()

    # ZVC
    point = "L4"
    plot_ZVC(3.0)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    ax.set_aspect("equal")
    plt.title("$J_{L4} < J < J_{L3}$")
    plt.grid("on",linestyle=":")


    #fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system", y=0.05)
    fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system")
    #plt.subplots_adjust(bottom=0.15,top=0.95)
    plt.savefig(f"{consts.PLOT_DIR}/CR3BP_jacobi_integral_variation.png")

def plot_anim():
    """
    builds animation of variation of jacobi energy and
    various ZVC variations
    """

    import matplotlib.animation as animation

    fig, ax = plt.subplots()

    x = np.linspace(consts.LOWERLIM, consts.UPPERLIM, 150)
    y = np.linspace(consts.LOWERLIM, consts.UPPERLIM, 150)

    X, Y = np.meshgrid(x,y)

    r1 = ((X+consts.MU)**2 + Y**2)**(1/2)
    r2 = ((X-(1-consts.MU))**2 + Y**2)**(1/2)

    J = np.linspace(s_utils.get_J_point("L1"),s_utils.get_J_point("L4"),150)

    # ZVC function
    def ZVC(i):
        return X**2 + Y**2 + (2*(1-consts.MU))/r1 + (2*consts.MU)/r2 - J[i]
    
    # initial contour
    global c
    c = ax.contour(X,Y,ZVC(0),[0],colors="black")
    #jlabel = ax.text(0.5,1.0, "", horizontalalignment="center",verticalalignment="center", transform=ax.transAxes)

    # animation function
    def animate(i):
        global c
        Z = ZVC(i)
        for col in c.collections:
            col.remove()
        c = ax.contour(X,Y,Z,[0],colors="black")
        ax.set_title(f"J = {round(J[i],4):.4f} for the Earth-Moon system")
        return c

    # libration points
    plot_libration_points()
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([consts.LOWERLIM,consts.UPPERLIM])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([consts.LOWERLIM,consts.UPPERLIM])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM+consts.STEP,step=consts.STEP))
    plt.ylabel("y (nondim)")
    ax.set_aspect("equal")
    plt.grid("on",linestyle=":")
    
    anim = animation.FuncAnimation(fig, animate, frames=len(J), repeat=False)
    anim.save(f"{consts.PLOT_DIR}/ZVC.mp4", writer="ffmpeg", fps=15)