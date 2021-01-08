"""
utilities for the main krill script, supporting static plotting
includes:
plotting zero velocity curve
plotting lagrange (libration) points
getting jacobi integral at a libration point


jacobi integral at L1 is 3.1885282305574663
jacobi integral at L2 is 3.1730187952481117
jacobi integral at L3 is 3.01215069525063
jacobi integral at L4 is 2.987993719716
L4 and L5 energy is same due to symmetry
J > J1 = 3.2
J2 < J < J1 = 3.180
J3 < J < J2 = 3.10
J4 < J < J3 = 2.99
"""

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import consts


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

def plot_lagrange_points():
    """
    computes the position of the 5 lagrange points
    https://www.mat.univie.ac.at/~westra/lagrangepoints.pdf
    """

    # L1 -- between the two masses
    # r1 = consts.MU, r2 = 1-consts.MU
    # x = a(1 - cbrt(M2/3M1)) -> x = (1-cbrt(consts.MU/3))
    L1 = [(1-math.pow(consts.MU/3,1/3)),0]
    #print(f"L1 at {L1}")
    plt.text(L1[0],L1[1],"L1")
    plt.scatter(L1[0],L1[1],s=5,c="black")

    # L2 -- on the line of the masses, beyond the smaller one
    # x = r2(1 + cbrt(M2/3M1)) -> x = (1-consts.MU)(1+cbrt(consts.MU/3)) 
    L2 = [(1-consts.MU)*(1+math.pow(consts.MU/3,1/3)),0]
    #print(f"L2 at {L2}")
    plt.text(L2[0],L2[1],"L2")
    plt.scatter(L2[0],L2[1],s=5,c="black")

    # L3 -- on the line of masses, beyond the larger one
    # x = r2(1 + 17/12 M2/M1) -> (1-consts.MU)(1 + 17/12 consts.MU)
    L3 = [-(1-consts.MU)*(1+17/12 * consts.MU),0]
    #print(f"L3 at {L3}")
    plt.text(L3[0],L3[1],"L3")
    plt.scatter(L3[0],L3[1],s=5,c="black")

    # L4 -- corner of top equilateral triangle
    # eq 2.45 
    # https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Dissertations/2011_CraigDavis.pdf
    # x = 1/2 - consts.MU, y = sqrt3/2
    L4 = [1/2-consts.MU,math.sqrt(3)/2]
    #print(f"L4 at {L4}")
    plt.text(L4[0],L4[1],"L4")
    plt.scatter(L4[0],L4[1],s=5,c="black")

    # L5 -- corner of bottom equilateral triangle
    L5 = [1/2-consts.MU,-math.sqrt(3)/2]
    #print(f"L5 at {L5}")
    plt.text(L5[0],L5[1],"L5")
    plt.scatter(L5[0],L5[1],s=5,c="black")

def get_J(point):
    """
    gets the jacobi integral at a specified lagrange point

    args:
        point -- string (L1-L5)
    
    returns:
        J -- jacobi integral
    """

    if point == "L1":
        rx = 1-math.pow(consts.MU/3,1/3)
        ry = 0
    elif point == "L2":
        rx = (1-consts.MU)*(1+math.pow(consts.MU/3,1/3))
        ry = 0
    elif point == "L3":
        rx = -(1-consts.MU)*(1+17/12 * consts.MU)
        ry = 0
    elif point == "L4":
        rx = 1/2-consts.MU
        ry = math.sqrt(3)/2 
    elif point == "L5":
        rx = 1/2-consts.MU
        ry = -math.sqrt(3)/2 
    else:
        print(f"no point exists for argument {point}")
        exit(-1)

    r1 = math.sqrt(math.pow((rx+consts.MU),2) + math.pow(ry,2))
    r2 = math.sqrt(math.pow((rx-(1-consts.MU)),2) + math.pow(ry,2))

    J = math.pow(rx,2) + math.pow(ry,2) + (2*(1-consts.MU))/r1 + (2*consts.MU)/r2

    return J