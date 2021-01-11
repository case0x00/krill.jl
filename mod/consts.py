"""
constants
"""

# plot limits
LOWERLIM = -1.2
UPPERLIM = 1.4
STEP = 0.2
PLOT_DIR = "plots"

def set_mu(mu):
    """
    set the mass parameter to be used in other files
    """
    print(f"SET MU AS {mu}")
    global MU
    MU = mu

def init_conds(rx0,ry0,vx0,vy0):
    """
    sets the initial conditions
    """

    global RX0,RY0,VX0,VY0
    RX0 = rx0
    RY0 = ry0
    VX0 = vx0
    VY0 = vy0