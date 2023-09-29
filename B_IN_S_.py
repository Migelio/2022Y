import numpy as np
import matplotlib.pyplot as plt
from BINS_class import Bins
fi0 = 55.811527 * np.pi / 180
la0 = 38.5 * np.pi / 180
h = 170
bins = Bins(fi0, la0, h)
file = open(r'C:\myProgram\py\kurs4\data.txt', 'r')
fly = file.readlines()
phi_errors = []
la_errors = []
f = []

for i in range(len(fly)):
    fly[i] = fly[i].split('\t')

    Omx, Omy, Omz, Ax, Ay, Az, pitch, roll, head, Vn, Ve = fly[i]
    
    bins.OrientStep(float(Omx)/3600/180*np.pi, float(Omy)/3600/180*np.pi, float(Omz)/3600/180*np.pi )
    
    phi, la, h, Vx, Vy, Vz = bins.NaviStep(np.array([float(Ax), float(Ay), float(Az)]))
    phi_errors.append((phi - fi0) * 6371000)
    la_errors.append((la - la0) * 6371000 * np.cos(phi))

    Rz = Bins.b / np.sqrt(1 - Bins.e**2 * np.cos(phi)**2)

    Omega_1 = np.array([[0, -Omz, Omy], 
                            [Omz, 0, -Omx], 
                            [-Omy, Omx, 0]])
        
    Vector_2 = np.array([Bins.u*np.cos(phi), Bins.u*np.sin(phi), 0]) + np.array([Vz/((Rz+h)), -Vz*np.sin(phi)/((Rz+h)*np.cos(phi)), 0]) + np.array([0, 0, -Vx/(Rz + h)])
        
    Omega_2 = np.array([[0, -Vector_2[2], Vector_2[1]], 
                            [Vector_2[2], 0, -Vector_2[0]], 
                            [-Vector_2[1], Vector_2[0], 0]])
    dC = C.dot(Omega_1) + Omega_2.dot(C)
    C = C + dC * Bins.dT

plt.plot(la_errors, phi_errors, color = 'black')
plt.grid()
plt.show()

plt.plot(range(len(fly)), phi_errors, color = 'black')
plt.grid()
plt.show()

plt.plot(range(len(fly)), la_errors, color = 'black')
plt.grid()
plt.show()
    
