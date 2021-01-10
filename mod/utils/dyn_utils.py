"""
utils for the main krill script, supporting animating the
third body

"""

import math
import tqdm

import mod.consts as consts

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
    #vx = state[2]
    #vy = state[3]

    # position of third body relative to the other two
    r1 = math.sqrt(math.pow((rx+consts.MU),2) + math.pow(ry,2))
    r2 = math.sqrt(math.pow((rx-(1-consts.MU)),2) + math.pow(ry,2))
    #print(f"r1 {r1}; r2 {r2}")

    # eq 2.18
    #ax = 2*vy + rx - ((1-consts.MU)*(rx+consts.MU))/(math.pow(r1,3)) - (consts.MU*(rx-(1-consts.MU)))/(math.pow(r2,3))
    #ay = -2*vx + y - ((1-consts.MU)*ry)/(math.pow(r1,3)) - (consts.MU*ry)/(math.pow(r2,3))

    ax = - ((1-consts.MU)*(rx+consts.MU))/(math.pow(r1,3)) - (consts.MU*(rx-(1-consts.MU)))/(math.pow(r2,3))
    ay = - ((1-consts.MU)*ry)/(math.pow(r1,3)) - (consts.MU*ry)/(math.pow(r2,3))

    #print(f"ax {ax}; ay {ay}")


    # jacobi integral -> do we need this?
    #J = math.pow(rx,2) + math.pow(ry,2) + (2*(1-consts.MU))/r1 + (2*consts.MU)/r2 - (math.pow(vx,2) + math.pow(vy,2))

    return [ax,ay]
