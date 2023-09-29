import numpy as np
from numpy import sin, cos

class Bins:
    ge = 9.78049
    a = 6378137.0
    b = 6356752.3142
    u = 0.7292115e-4
    e = np.sqrt(a**2 + b**2) / a
    dT = 1
    detC = 1
    def __init__(self, fi0, la0, h):
        self.phi = fi0
        self.la = la0
        self.h = h

        self.data = []
        self.file = open(r'C:\myProgram\py\kurs4\data.txt', 'r')
        self.fly = self.file.readlines()
        for i in range(len(self.fly)):
            self.fly[i] = self.fly[i].split('\t')
        self.Ac = []
        self.Omc = []
        for i in range(1000):
            self.Omc.append([float(self.fly[i][0]) / 180 * np.pi / 3600, float(self.fly[i][1]) / 180 * np.pi / 3600, float(self.fly[i][2]) / 180 * np.pi / 3600])
            self.Ac.append([float(self.fly[i][3]), float(self.fly[i][4]), float(self.fly[i][5])])
        self.Omc = np.mean(np.array(self.Omc), axis = 0)
        self.Ac = np.mean(np.array(self.Ac), axis = 0)
        self.Nc = np.cross(self.Ac, self.Omc)
        self.Matrix_c = np.array([self.Ac, self.Omc, self.Nc]).transpose()
        self.An = self.Gravi(fi0, h)
        self.Omn = np.array([Bins.u * np.cos(fi0), Bins.u * np.sin(fi0), 0])
        self.Nn = np.cross(self.An, self.Omn)
        self.Matrix_n = np.array([self.An, self.Omn, self.Nn]).transpose()
        self.C = self.Matrix_n.dot(np.linalg.inv(self.Matrix_c))

        self.Vx = 0
        self.Vy = 0
        self.Vz = 0
        self.VC = []
        for i in range(3):
            self.VC.append(np.zeros((3,3)))

    def Gravi(self, fi0, h):
        q = Bins.u ** 2 * Bins.a / Bins.ge
        e = ((Bins.a ** 2 - Bins.b ** 2) / Bins.a ** 2) ** 0.5
        alpha = (Bins.a - Bins.b) / Bins.a
        betta  = 5/2 * q - alpha - 17/14 * q * alpha
        betta1 = alpha  ** 2 / 8 - 5 / 8 * q * alpha
        g0 = Bins.ge * (1 + betta * sin(fi0) ** 2 + betta1 * sin(fi0) ** 2)
        gz = 0
        gx = 0.5 * Bins.ge * h / Bins.a * (4 * q - e ** 2) * sin(2 * fi0)
        gy = g0 + h / Bins.a * (3 * h / Bins.a - 2*q * cos(fi0) **2 + e ** 2 * ( 3* sin(fi0) ** 2 - 1) - q * (1 + 6 * sin(fi0) ** 2))
        return np.array([gx, gy, gz])
    
    def OrientStep(self, Omx, Omy, Omz):
        
        Rz = Bins.b / np.sqrt(1 - Bins.e**2 * np.cos(self.phi)**2)

        Omega_1 = np.array([[0, -Omz, Omy], 
                            [Omz, 0, -Omx], 
                            [-Omy, Omx, 0]])
        
        Vector_2 = np.array([Bins.u*np.cos(self.phi), Bins.u*np.sin(self.phi), 0]) + np.array([self.Vz/((Rz+self.h)), -self.Vz*np.sin(self.phi)/((Rz+self.h)*np.cos(self.phi)), 0]) + np.array([0, 0, -self.Vx/(Rz + self.h)])
        
        Omega_2 = np.array([[0, -Vector_2[2], Vector_2[1]], 
                            [Vector_2[2], 0, -Vector_2[0]], 
                            [-Vector_2[1], Vector_2[0], 0]])    
                       
        dC = self.C.dot(Omega_1) + Omega_2.dot(self.C)
        self.C = self.C + dC * Bins.dT
        # Принудительная ортогонализация
        E = 0.5 * (self.C.dot(np.transpose(self.C)) - np.ones((3, 3)))
        self.C = E.dot(self.C)
        DetC = self.C
        return self.C
    
    
                
    def NaviStep(self, A):
        Ax, Ay, Az = self.C.dot(A)
        Rz = Bins.b / np.sqrt(1 - Bins.e**2 * np.cos(self.phi)**2)
        g = self.Gravi(self.phi, self.h)
        phi_dot = -self.Vx / (Rz + self.h)
        la_dot = self.Vz / ((Rz + self.h)*np.cos(self.phi))
        Vx_dot = -(2*Bins.u + la_dot)*np.sin(self.phi)*self.Vz - phi_dot*self.Vy + g[0] + Ax
        Vy_dot = -(2*Bins.u + la_dot)*np.cos(self.phi)*self.Vz - phi_dot*self.Vx + g[1] + Ay
        Vz_dot = -(2*Bins.u*np.cos(self.phi) + la_dot*np.sin(self.phi))*self.Vy + (2*Bins.u + la_dot)*np.sin(self.phi)*self.Vx + g[2] + Az
        self.Vx += Vx_dot * Bins.dT
        self.Vy += Vy_dot * Bins.dT
        self.Vz += Vz_dot * Bins.dT
        phi_dot = -self.Vx / (Rz + self.h)
        la_dot = self.Vz / ((Rz + self.h)*np.cos(self.phi))
        self.phi += phi_dot * Bins.dT
        self.la += la_dot * Bins.dT
        return (self.phi, self.la, self.h, self.Vx, self.Vy, self.Vz)
