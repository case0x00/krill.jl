#!/usr/bin/env python3

"""
simulates solutions to the CR3BP, allowing for variable
mass parameter, initial conditions, and jacobi integral
in the synodic reference frame (barycenter is at the origin)
thus angular velocity n = 1

m1 at (-mu, 0) and m2 at (1-mu, 0) so the distance is unity
orbital period is 2pi
m3 is r1 and r2 away from m1 and m2 respectively

times units of ~104h (one sidereal month /2)
length units of 384,400 km
velocity units of 1024 m/s

the hamiltonian (total energy) is
H_{CR3BP} = 1/2 (p_x^2 + p_y^2) + yp_x - xp_y - (1-mu)/r1 - (mu)/r2
where px = xdot-y, py = ydot+x, r1^2=(x-xe)^2+y^2, r2^2 = (x-xm)^2+y^2
H_{CR3BP} = -C/2 in the synodic system

args:
    mu -- mass parameter
    rx0 -- initial x
    ry0 -- initial y
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

# mass parameter
MU = float(sys.argv[1])

def f(state):
    """
    function to apply RK4 to, which updates the acceleration
    vectors for the third body orbiting the other two

    args:
        state -- [rx,ry,vx,vy]
    returns:
        state -- [ax,ay]

    """

    # body positions on x-axis relative to mass ratios
    body1_rx = -MU
    body2_rx = 1-MU

    # positions of third body
    rx = state[0]
    ry = state[1]

    # positions relative to first and second body
    r1 = math.sqrt((rx-body1_rx) ** 2 + ry ** 2)
    r2 = math.sqrt((rx-body2_rx) ** 2 + ry ** 2)

    # velocity of the third body
    #vx = state[2]
    #vy = state[3]

    # accelerations
    #ax = 2 * vy + rx -(1-MU) * (rx-body1_rx)/(r1 ** 3) - MU * (rx-body2_rx)/(r2 ** 3)
    #ay = -2 * vy + ry -(1-MU) * ry/(r1 ** 3) - MU * ry/(r2 ** 3)

    ax = -(1-MU) * (rx-body1_rx)/(r1 ** 3) - MU * (rx-body2_rx)/(r2 ** 3)
    ay = -(1-MU) * ry/(r1 ** 3) - MU * ry/(r2 ** 3)

    return [ax, ay]


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



def plot_lagrange_points():
    """
    computes the position of the 5 lagrange points
    https://www.mat.univie.ac.at/~westra/lagrangepoints.pdf
    idk how useful it is given the assumptions it makes

    bearing in mind, M1 (-mu,0) > M2 (1-mu,0)
    and the mass param mu is M2/M1

    returns:
        [L1(x,y),L2(x,y),L3(x,y),L4(x,y),L5(x,y)]
    """

    # L1 -- between the two masses
    # r1 = mu, r2 = 1-mu
    # x = a(1 - cbrt(M2/3M1)) -> x = (1-cbrt(MU/3))
    L1 = [(1-math.pow(MU/3,1/3)),0]
    print(f"L1 at {L1}")
    plt.text(L1[0],L1[1],"L1")

    # L2 -- on the line of the masses, beyond the smaller one
    # x = r2(1 + cbrt(M2/3M1)) -> x = (1-mu)(1+cbrt(mu/3)) 
    L2 = [(1-MU)*(1+math.pow(MU/3,1/3)),0]
    print(f"L2 at {L2}")
    plt.text(L2[0],L2[1],"L2")

    # L3 -- on the line of masses, beyond the larger one
    # x = r2(1 + 17/12 M2/M1) -> (1-mu)(1 + 17/12 mu)
    L3 = [-(1-MU)*(1+17/12 * MU),0]
    print(f"L3 at {L3}")
    plt.text(L3[0],L3[1],"L3")

    # L4 -- corner of top equilateral triangle
    # dist between the two bodies is 1, so each triangle side is 0.5?
    # triangle is 1^2 = 0.5^2 + h^2 -> h = sqrt(1^2-0.5^2)
    h = math.sqrt(1**2 - 0.5**2)
    L4 = [0.5,h]
    print(f"L4 at {L4}")
    plt.text(L4[0],L4[1],"L4")



    # L5 -- corner of bottom equilateral triangle
    L5 = [0.5,-h]
    print(f"L5 at {L5}")
    plt.text(L5[0],L5[1],"L5")





#Def RK4(state):
#    """
#    uses the runge-kutta 4 method to approximate the solutions
#    to the equations of motion for the 3rd body and returns the 
#    x and y positions.


#    applying RK4 to both equations for acceleration at the same time (NOT JOINT)
#    it increments each respective acceleration

#    args:
#        state -- [rx,ry,vx,vy]
#    returns:
#        state -- [rx,ry]
#    """

#    # stop condition as duration [s]

#    duration = 2 * 365 * 24 * 3600
#    # number of steps
#    n_steps = 10000
#    # step size
#    h = int(duration/n_steps)

#    # time range
#    t = range(0, duration, h)
#    t_i = 0
#    t_max = max(t)

#    rx = []
#    ry = []


#    i = 0
#    stop = False
#    while stop == False:
#        if t_i == t_max or t_i > t_max:
#            stop = True
#        elif i > n_steps:
#            stop = True
#        else:
#            # RK4
#            if i%100 == 0:
#                print(f"it {i}, t_i {t_i}, h {h}")
#            rx.append(state[0])
#            ry.append(state[1])


#            k1 = h * get_state(state)
##            k2 = h * get_state(state + 1/2 * k1)
#            k2 = h * get_state(state + [1/2 * x for x in k1])
##            k3 = h * get_state(state + 1/2 * k2)
#            k3 = h * get_state(state + [1/2 * x for x in k2])
#            k4 = h * get_state(state + k3)
#            state = [1/6 * x for x in (k1 + 2 * k2 + 2 * k3 + k4)]
#            print(len(state))

#            t_i += h
#            i += 1

#    print("end loop")
#    return [rx,ry]


#def animate(i):
#    """
#    animate function for matplotlib.FuncAnimation
#    """
#
#    data = 


def main():
    x0 = float(sys.argv[2])
    y0 = float(sys.argv[3])
    vx0 = float(sys.argv[4])
    vy0 = float(sys.argv[5])

#    tic = time.perf_counter()
#    x,y = RK4([x0,y0,vx0,vy0])
#    toc = time.perf_counter()
#    print(f"completed in {toc-tic:0.4f}s")

#    t = len(x)


    # plotting
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.html
    # https://matplotlib.org/3.1.0/gallery/color/named_colors.html

    matplotlib.rcParams['font.size'] = 10.
    matplotlib.rcParams['font.family'] = 'Serif'
    matplotlib.rcParams['axes.labelsize'] = 10.
    matplotlib.rcParams['xtick.labelsize'] = 10.
    matplotlib.rcParams['ytick.labelsize'] = 10.


    fig, ax = plt.subplots(figsize=(12,10))

    # lagrange points
    plot_lagrange_points()

    # body 3
#    plt.scatter(x,y, marker="*")
    # body 1 (larger)
    plt.scatter(-MU,0,c="royalblue",s=60)
    print(f"body 1 at {[-MU,0]}")
    # body 2 (smaller)
    plt.scatter(1-MU,0,c="slategray",s=20)
    print(f"body 2 at {[1-MU,0]}")
    plt.xlim([-1,1])
    plt.xticks(np.arange(-1.2,1.4,step=0.2))
    plt.ylim([-1,1])
    plt.yticks(np.arange(-1.2,1.4,step=0.2))
    plt.grid("on",linestyle=":")
    plt.title(f"CR3BP in the synodic frame in the earth-moon system")
    plt.savefig("CR3BP.png")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("incorrect number of arguments") 
        exit(-1)

    main()