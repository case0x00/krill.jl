#!/usr/bin/env python3

"""
simulates solutions to the CR3BP, allowing for variable
mass parameter, initial conditions, and jacobi integral
in the synodic reference frame (barycenter is at the origin)
thus angular velocity n = 1

args:
    mu -- mass parameter
    rx0 -- initial rx (of third body)
    ry0 -- initial ry
    vx0 -- initial vx
    vy0 -- initial vy

"""

import sys
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
from tqdm import tqdm

import consts
import static_utils as static
import dyn_utils as dyn

# mass parameter of m2/m* = m2/(m1+m2)
consts.set_mu(float(sys.argv[1]))

PROGRAM = "static_anim"


def RK4(initial_state):
    """
    uses the runge-kutta 4 method to approximate the solutions to
    the equations of motion for the third body and returns the x,y
    positions
    https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

    args:
        initial_state -- [rx0,ry0,vx0,vy0]
    """

    #nsteps = 10000
    #tstep = 10

    duration = 10 * 365 * 24 * 3600
    nsteps = 10000
    tstep = int(duration/nsteps)

    # these should always be linear values
    rx = initial_state[0]
    ry = initial_state[1]
    vx = initial_state[2]
    vy = initial_state[3]
    ax = 0
    ay = 0
    rx_out = [rx]
    ry_out = [ry]

    for i in tqdm(range(0, duration, nsteps)):
        # iteratively apply RK4

        # calling f(state) returns [ax,ay]

        # updating acceleration
        # k1 = f(rx,ry,vx,vy) -> [ax,ay]
        k1 = f([rx,ry])

        # k2 = f(x+h/2, y+h/2*k1) -> [ax,ay]
        # apply multipler to rx, ry

        k2 = f([rx+tstep/2*k1[0],ry+tstep/2*k1[1]])

        # k3 = f(x+h/2, y+h/2*k2) -> [ax,ay]
        k3 = f([rx+tstep/2*k2[0],ry+tstep/2*k2[1]])

        # k4 = f(x+h, y+h*k3) -> [ax,ay]
        k4 = f([rx+tstep*k3[0],ry+tstep*k3[1]])

        # yn+1 = yn + 1/6 * h (k1 + 2k2 + 2k3 + k4)
        # updates ax, ay
        ax += 1/6 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        ay += 1/6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])

        # update velocity
        vx += ax * tstep
        vy += ay * tstep

        # update position
        rx += vx * tstep
        ry += vy * tstep
        rx_out.append(rx)
        ry_out.append(ry)


    return [rx_out,ry_out]


def f(state):
    """
    the CR3BP equations of motion. this is the function to
    apply RK4 to, which updates the acceleration
    vectors for the third body orbiting the other two

    args:
        state -- [rx,ry,vx,vy]
    returns:
        state -- [ax,ay]
    """

    rx = state[0]
    ry = state[1]
    vx = state[2]
    vy = state[3]

    # position of third body relative to the other two
    r1 = math.sqrt(math.pow((rx+consts.MU),2) + math.pow(ry,2))
    r2 = math.sqrt(math.pow((rx-(1-consts.MU)),2) + math.pow(ry,2))
    print(f"r1 {r1}; r2 {r2}")

    # eq 2.18
    ax = 2*vy + rx - ((1-consts.MU)*(rx+consts.MU))/(math.pow(r1,3)) - (consts.MU*(rx-(1-consts.MU)))/(math.pow(r2,3))
    ay = -2*vx + y - ((1-consts.MU)*ry)/(math.pow(r1,3)) - (consts.MU*ry)/(math.pow(r2,3))

    # jacobi integral -> do we need this?
    J = math.pow(rx,2) + math.pow(ry,2) + (2*(1-consts.MU))/r1 + (2*consts.MU)/r2 - (math.pow(vx,2) + math.pow(vy,2))

    return [ax,ay]


def plot_synodic_gif(J):
    #x0 = float(sys.argv[2])
    #y0 = float(sys.argv[3])
    #vx0 = float(sys.argv[4])
    #vy0 = float(sys.argv[5])

    #tic = time.perf_counter()
    #x,y = RK4([x0,y0,vx0,vy0])
    #toc = time.perf_counter()
    #print(f"completed in {toc-tic:0.4f}s")

    #t = len(x)

    # lagrange points
    #plot_lagrange_points()
    #plot_ZVC(J)





    pass



def plot_static():
    """
    builds 4 subplots with variations in jacobi integral
    to demonstrate the ZVC
    """

    fig = plt.figure(figsize=(12,10))

    ax = fig.add_subplot(221)

    # lagrange points
    static.plot_lagrange_points()

    # ZVC
    point = "L1"
    static.plot_ZVC(3.2)


    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J > J_{L1}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(222)

    # lagrange points
    static.plot_lagrange_points()

    # ZVC
    point = "L2"
    static.plot_ZVC(3.180)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L2} < J < J_{L1}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(223)

    # lagrange points
    static.plot_lagrange_points()

    # ZVC
    point = "L3"
    static.plot_ZVC(3.10)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L3} < J < J_{L2}$")
    plt.grid("on",linestyle=":")


    ax = fig.add_subplot(224)

    #lagrange points
    static.plot_lagrange_points()

    # ZVC
    point = "L4"
    static.plot_ZVC(3.0)
 
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L4} < J < J_{L3}$")
    plt.grid("on",linestyle=":")


    #fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system", y=0.05)
    fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system")
    #plt.subplots_adjust(bottom=0.15,top=0.95)
    plt.savefig("plots/CR3BP_jacobi_integral_variation.png")




def plot_static_anim():
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

    J = np.linspace(static.get_J("L1"),static.get_J("L4"),150)

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
        ax.set_title(f"Jacobi integral J = {round(J[i],4):.4f} for the Earth-Moon system")
        return c

    # lagrange points
    static.plot_lagrange_points()
    # body 1 (larger)
    plt.scatter(-consts.MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-consts.MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(consts.LOWERLIM,consts.UPPERLIM,step=consts.STEP))
    plt.ylabel("y (nondim)")
    plt.grid("on",linestyle=":")
    
    anim = animation.FuncAnimation(fig, animate, frames=len(J), repeat=False)
    anim.save("anims/ZVC.mp4", writer="ffmpeg", fps=15)



def main():
    # plotting
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.html
    # https://matplotlib.org/3.1.0/gallery/color/named_colors.html

    matplotlib.rcParams['font.size'] = 10.
    matplotlib.rcParams['font.family'] = 'Serif'
    matplotlib.rcParams['axes.labelsize'] = 10.
    matplotlib.rcParams['xtick.labelsize'] = 10.
    matplotlib.rcParams['ytick.labelsize'] = 10.

    tic = time.perf_counter()

    if PROGRAM == "static":
        print("BUILDING STATIC PLOT OF ZVC")
        plot_static()
    elif PROGRAM == "static_anim":
        print("BUILDING ANIMATION OF ZVC")
        plot_static_anim()
    elif PROGRAM == "dynamic_anim":
        print("NOT YET IMPLEMENTED")
#    plot_synodic_gif(static.get_J("L2", consts.MU))
    else:
        print("PROGRAM not specified")
        exit(-1)

    toc = time.perf_counter()
    print(f"COMPLETED IN {toc-tic:0.4f}s")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("INCORRECT NUMBER OF ARGUMENTS") 
        exit(-1)

    main()