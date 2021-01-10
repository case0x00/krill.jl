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
    global MU
    MU = mu

def init_conds(rx0,ry0,xv0,xy0):
    """
    sets the initial conditions
    """

    global RX0,RY0,XV0,XY0
    RX0 = rx0
    RY0 = ry0
    XV0 = xv0
    XY0 = xy0