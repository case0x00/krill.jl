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
import os
import time

current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

import mod.consts as consts
import mod.plot.params as p
import mod.plot.static as s
import mod.plot.dynamic as d


# mass parameter of m2/m* = m2/(m1+m2)
consts.set_mu(float(sys.argv[1]))

# set initial conditions
consts.init_conds(float(sys.argv[2]),
                  float(sys.argv[3]),
                  float(sys.argv[4]),
                  float(sys.argv[5]))


# program to execute
PROGRAM = "dynamic_anim"

def main():
    # initialize matplotlib params
    p.init_params()

    tic = time.perf_counter()

    if PROGRAM == "static":
        print("BUILDING STATIC PLOT OF ZVC")
        s.plot()
    elif PROGRAM == "static_anim":
        print("BUILDING ANIMATION OF ZVC")
        s.plot_anim()
    elif PROGRAM == "dynamic_anim":
        print("BUILDING ANIMATION OF THIRD BODY")
        d.plot_anim()
    else:
        print("PROGRAM not specified")
        exit(-1)

    toc = time.perf_counter()
    print(f"COMPLETED IN {toc-tic:.4f}s")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("INCORRECT NUMBER OF ARGUMENTS") 
        exit(-1)

    main()