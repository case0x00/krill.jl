"""
utils for the main krill script, supporting animating the
third body

"""

import math
import tqdm
from numpy import linalg

import mod.consts as consts

def old_RK4(initial_state):
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


def old_f(state):
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


def RK4(initial_state):
    """

    """

    rx = initial_state[0]
    ry = initial_state[1]
    vx = initial_state[2]
    vy = initial_state[3]

    x = [rx,ry]

    # each k value represents
    # [vx,vy,ax,ay]


    # each f call requires a state, and returns the
    # estimated state derivative
    # h is a timestep, and by multiplying the
    # state derivative by the timestep, it basically
    # integrates and you add the old state after 
    # 1 timestep and the new state
    
    # k1 = f(t, s(rx,ry,vx,vy))
    k1 = f(t,x)

    # k2 = f(t+dt, s(rx,ry,vx,vy) + dt*old_s(vx,vy,ax,ay))
    k2 = f(t+h/2, x+h/2*k1)

    k3 = f(t+h/2, x+h/2*k2)

    k4 = f(t+h, x+h*k3)

    # xf = s(rx,ry,vx,vy) + 1/6*(old_s(vx,vy,ax,ay))
    # xf(rx,ry,vx,vy)
    xf = x + 1/6*(k1 + 2*k2 + 2*k3 + k4)

    # for each timestep, the RK4 solves for xf
    # xf value encodes the updated position and velocites, which are fed back into the RK4
    # to update the next position






def f(t,state):
    """
    the CR3BP equations of motion

    args:
        state -- [rx,ry,vx,vy]
    returns:
        state -- [vx,vy,ax,ay]


    i need a time variable

    """
    a = [0,0]

    r = [state[0],state[1]]
    v = [state[2],state[3]]

    rb1 = [-consts.MU,0]
    rb2 = [1-consts.MU,0]

    r1 = [a-b for a,b in zip(r,rb1)]
    r2 = [a-b for a,b in zip(r,rb2)]

    r13 = linalg.norm(r1,2)**3
    r32 = linalg.norm(r2,2)**3

    c1 = (1-consts.MU)/r13
    c2 = consts.MU/r23

    a[0] = 2 * v[1] + r[0] - c1 * (r[0] + consts.MU) - c2 * (r[0] - (1-consts.MU))
    a[1] = -2 * v[0] + r[1] - c1 * r[1] - c2 * r[1]

    # position -- r(1x2) as input
    # velocity -- v(1x2) as input
    # acc -- a(1x2) as output

    return [r[0],r[1],a[0],a[1]]



