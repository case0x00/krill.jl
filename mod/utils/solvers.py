"""
solving differential equations
"""

def RK4(f, rx, ry, vx, vy, h, i):
    """
    applies the RK4 to approximate solutions to the
    equations of motion
    """

    k1 = f(rx,ry,vx,vy)
    k2 = f(rx + (h * k1[0])/2, ry + (h * k1[1])/2, vx + (h * k1[2])/2, vy + (h * k1[3])/2)
    k3 = f(rx + (h * k2[0])/2, ry + (h * k2[1])/2, vx + (h * k2[2])/2, vy + (h * k2[3])/2)
    k4 = f(rx + h * k3[0], ry + h * k3[1], vx + h * k3[2], vy + h * k3[3])

    if i < 10:
        print(f"i {i}")
        print(f"k1 {k1} -> rx {rx}; ry {ry}")
        print(f"k2 {k2} -> rx {rx}; ry {ry}")
        print(f"k3 {k3} -> rx {rx}; ry {ry}")
        print(f"k4 {k4} -> rx {rx}; ry {ry}")


    rx_step = rx + (h/6) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])    
    ry_step = ry + (h/6) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    vx_step = vx + (h/6) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    vy_step = vy + (h/6) * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])

    return rx_step, ry_step, vx_step, vy_step