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

    fi = (90 - 180 / 22 * i)*pi/180
    lambd = (360 - 360 / 22 * i)*pi/180
    H = 3000 / 22 * i
    print('fi lambd H', fi, lambd, H)

    yaw = (360 - 360/22 * i)*pi/180  # psi
    pitch = (90 - 180/22 * i)*pi/180  # ve
    roll = (180 - 360/22 * i)*pi/180  # gamm
    print('ypr', yaw, pitch, roll)
    H = 3000 / 22 * i

    An = array([[cos(yaw)*cos(pitch),  sin(yaw)*sin(roll)-cos(yaw)*sin(pitch)*cos(roll),   sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)*sin(roll)],
                [sin(pitch),           cos(
                    pitch)*cos(roll),                              -cos(pitch)*sin(roll)],
                [-sin(yaw)*cos(pitch), cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll),   cos(yaw)*cos(roll)-sin(yaw)*sin(pitch)*sin(roll)]])
    print('\ndet(An)\n', np.linalg.det(An), '\nAn\n', An)
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
        d_lambd = Vg[2] / (Rz + H) * cos(fi) * dt
        Ht = H + Vg[1] * dt
        fi_t = (fi + d_fi)
        lambd_t = (lambd + d_lambd)
    print('H', fi, lambd, H)
    print('Ht', fi_t, lambd_t, Ht)
    print((fi_t-fi)*Rz*H)
    print((lambd_t-lambd))
    print(Ht-H)


dano(8)
