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
from matplotlib.animation import FuncAnimation
import numpy as np
import time
from tqdm import tqdm

# mass parameter of m2/m* = m2/(m1+m2)
MU = float(sys.argv[1])

# plot limits
LOWERLIM = -1.2
UPPERLIM = 1.4
STEP = 0.2

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
    r1 = math.sqrt(math.pow((rx+MU),2) + math.pow(ry,2))
    r2 = math.sqrt(math.pow((rx-(1-MU)),2) + math.pow(ry,2))
    print(f"r1 {r1}; r2 {r2}")

    # eq 2.18
    ax = 2*vy + rx - ((1-MU)*(rx+MU))/(math.pow(r1,3)) - (MU*(rx-(1-MU)))/(math.pow(r2,3))
    ay = -2*vx + y - ((1-MU)*ry)/(math.pow(r1,3)) - (MU*ry)/(math.pow(r2,3))

    # jacobi integral -> do we need this?
    J = math.pow(rx,2) + math.pow(ry,2) + (2*(1-MU))/r1 + (2*MU)/r2 - (math.pow(vx,2) + math.pow(vy,2))

    return [ax,ay]


def plot_lagrange_points():
    """
    computes the position of the 5 lagrange points
    https://www.mat.univie.ac.at/~westra/lagrangepoints.pdf
    """

    # L1 -- between the two masses
    # r1 = mu, r2 = 1-mu
    # x = a(1 - cbrt(M2/3M1)) -> x = (1-cbrt(MU/3))
    L1 = [(1-math.pow(MU/3,1/3)),0]
    #print(f"L1 at {L1}")
    plt.text(L1[0],L1[1],"L1")
    plt.scatter(L1[0],L1[1],s=5,c="black")

    # L2 -- on the line of the masses, beyond the smaller one
    # x = r2(1 + cbrt(M2/3M1)) -> x = (1-mu)(1+cbrt(mu/3)) 
    L2 = [(1-MU)*(1+math.pow(MU/3,1/3)),0]
    #print(f"L2 at {L2}")
    plt.text(L2[0],L2[1],"L2")
    plt.scatter(L2[0],L2[1],s=5,c="black")

    # L3 -- on the line of masses, beyond the larger one
    # x = r2(1 + 17/12 M2/M1) -> (1-mu)(1 + 17/12 mu)
    L3 = [-(1-MU)*(1+17/12 * MU),0]
    #print(f"L3 at {L3}")
    plt.text(L3[0],L3[1],"L3")
    plt.scatter(L3[0],L3[1],s=5,c="black")

    # L4 -- corner of top equilateral triangle
    # eq 2.45 
    # https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf
    # x = 1/2 - mu, y = sqrt3/2
    L4 = [1/2-MU,math.sqrt(3)/2]
    #print(f"L4 at {L4}")
    plt.text(L4[0],L4[1],"L4")
    plt.scatter(L4[0],L4[1],s=5,c="black")

    # L5 -- corner of bottom equilateral triangle
    L5 = [1/2-MU,-math.sqrt(3)/2]
    #print(f"L5 at {L5}")
    plt.text(L5[0],L5[1],"L5")
    plt.scatter(L5[0],L5[1],s=5,c="black")


def plot_ZVC(J):
    """
    plots the hill region (zero velocity curves).

    https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf
    eq 2.53

    args:
        J -- jacobi integral

    """

    x = np.linspace(LOWERLIM, UPPERLIM, 150)
    y = np.linspace(LOWERLIM, UPPERLIM, 150)

    X, Y = np.meshgrid(x,y)

    r1 = ((X+MU)**2 + Y**2)**(1/2)
    r2 = ((X-(1-MU))**2 + Y**2)**(1/2)

    Z = X**2 + Y**2 + (2*(1-MU))/r1 + (2*MU)/r2 - J

    plt.contour(X,Y,Z,[0],colors="black")


def get_j(point):
    """
    gets the jacobi integral at a specified lagrange point

    args:
        point -- string (L1-L5)
    
    returns:
        J -- jacobi integral
    """

    if point == "L1":
        rx = 1-math.pow(MU/3,1/3)
        ry = 0
    elif point == "L2":
        rx = (1-MU)*(1+math.pow(MU/3,1/3))
        ry = 0
    elif point == "L3":
        rx = -(1-MU)*(1+17/12 * MU)
        ry = 0
    elif point == "L4":
        rx = 1/2-MU
        ry = math.sqrt(3)/2 
    elif point == "L5":
        rx = 1/2-MU
        ry = -math.sqrt(3)/2 
    else:
        print(f"no point exists for argument {point}")
        exit(-1)

    r1 = math.sqrt(math.pow((rx+MU),2) + math.pow(ry,2))
    r2 = math.sqrt(math.pow((rx-(1-MU)),2) + math.pow(ry,2))

    J = math.pow(rx,2) + math.pow(ry,2) + (2*(1-MU))/r1 + (2*MU)/r2

    return J



def plot_synodic_gif():
    #x0 = float(sys.argv[2])
    #y0 = float(sys.argv[3])
    #vx0 = float(sys.argv[4])
    #vy0 = float(sys.argv[5])

    #tic = time.perf_counter()
    #x,y = RK4([x0,y0,vx0,vy0])
    #toc = time.perf_counter()
    #print(f"completed in {toc-tic:0.4f}s")

    #t = len(x)

    pass






def plot_j_variation():
    fig = plt.figure(figsize=(12,10))

    ax = fig.add_subplot(221)

    # lagrange points
    plot_lagrange_points()

    # ZVC
    point = "L1"
    plot_ZVC(3.2)
    # jacobi integral at L1 is 3.1885282305574663
    # jacobi integral at L2 is 3.1730187952481117
    # jacobi integral at L3 is 3.01215069525063
    # jacobi integral at L4 is 2.987993719716
    # L4 and L5 energy is same due to symmetry
    # J > J1 = 3.2
    # J2 < J < J1 = 3.180
    # J3 < J < J2 = 3.10
    # J4 < J < J3 = 2.99

    # body 1 (larger)
    plt.scatter(-MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J > J_{L1}$")
    plt.grid("on",linestyle=":")





    ax = fig.add_subplot(222)

    # lagrange points
    plot_lagrange_points()

    # ZVC
    point = "L2"
    plot_ZVC(3.180)
 
    # body 1 (larger)
    plt.scatter(-MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L2} < J < J_{L1}$")
    plt.grid("on",linestyle=":")




    ax = fig.add_subplot(223)

    # lagrange points
    plot_lagrange_points()

    # ZVC
    point = "L3"
    plot_ZVC(3.10)
 
    # body 1 (larger)
    plt.scatter(-MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L3} < J < J_{L2}$")
    plt.grid("on",linestyle=":")





    ax = fig.add_subplot(224)

    #lagrange points
    plot_lagrange_points()

    # ZVC
    point = "L4"
    plot_ZVC(3.0)
 
    # body 1 (larger)
    plt.scatter(-MU,0,c="royalblue",s=60)
    # body 2 (smaller)
    plt.scatter(1-MU,0,c="slategray",s=20)

    plt.xlim([-1,1])
    plt.xticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.xlabel("x (nondim)")
    plt.ylim([-1,1])
    plt.yticks(np.arange(LOWERLIM,UPPERLIM,step=STEP))
    plt.ylabel("y (nondim)")
    plt.title("$J_{L4} < J < J_{L3}$")
    plt.grid("on",linestyle=":")


    #fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system", y=0.05)
    fig.suptitle(f"Zero velocity curves for four values of $J$ in the Earth-Moon system")
    #plt.subplots_adjust(bottom=0.15,top=0.95)
    plt.savefig("CR3BP_jacobi_integral_variation.png")





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
    plot_j_variation()
    toc = time.perf_counter()
    print(f"completed in {toc-tic:0.4f}s")



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("incorrect number of arguments") 
        exit(-1)

    main()