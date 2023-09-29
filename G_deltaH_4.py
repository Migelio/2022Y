from math import pi, cos, sin, sqrt, atan
import numpy as np
from numpy import array, dot, transpose, empty
Rz = 6371000
time = 100


def dano(i):
    fi = (90 - 180 / 22 * i)
    lambd = (360 - 360 / 22 * i)
    H = 3000 / 22 * i
    print('fi lambd H', fi, lambd, H)

    Vxc = 0.5 + 3/22 * i
    Vyc = -3 + 4/22 * i
    Vzc = 4 - 3/22 * i
    print('Vc_all', Vxc, Vyc, Vzc)

    fi0 = (90 - 180 / 22 * i)*pi/180
    lambd0 = (360 - 360 / 22 * i)*pi/180
    H0 = 3000 / 22 * i
    fi = fi0
    lambd = lambd0
    H = H0
    print('fi lambd H', fi, lambd, H)

    yaw = (360 - 360/22 * i)*pi/180  # psi
    pitch = (90 - 180/22 * i)*pi/180  # ve
    roll = (180 - 360/22 * i)*pi/180  # gamm
    print('ypr', yaw, pitch, roll)
    H = 3000 / 22 * i
    print("h", H)
    An = array([[cos(yaw)*cos(pitch),  sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll),   sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)*sin(roll)],
                [sin(pitch),           cos(
                    pitch)*cos(roll),                              -cos(pitch)*sin(roll)],
                [-sin(yaw)*cos(pitch), cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll),   cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll)]])

    '''1) Vc = Acn*VC'''
    Vc = array([[Vxc], [Vyc], [Vzc]])
    Vn = An.dot(Vc)
    print('Vn\n', Vn)
    Rn = Vn.dot(time)
    print('Rn\n', Rn)
    Vg = array([[0.41018], [-3.58367], [-0.60932]])
    dt = 0
    for x in range(0, time, 1):
        dt += 1
        d_fi = Vg[0] / (Rz + H) * dt
        d_lambd = Vg[2] / ((Rz + H) * cos(fi)) * dt
        H = H + Vg[1]*1
        print("h", H)
        fi = (fi + d_fi)
        lambd = (lambd + d_lambd)
    print('H', fi, lambd, H)
    print('fi l H', fi0, lambd0, H0)
    print('fi-fi0', (fi-fi0)*pi*(Rz+H)/180)
    print('lambd-lambd0', (lambd-lambd0)*pi*(Rz+H)*cos(fi)/180)
    print(H-H0)


dano(8)
