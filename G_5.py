from matplotlib import pyplot
from math import pi
from numpy import sin, cos, array, arange

def RAD(grad):
	return grad / 180 * pi

def GRAD(grad):
	return grad * 180 / pi

def matr (psi, theta, gamma):
	Cy = array([[1,           0,          0],
	            [0,  cos(gamma), sin(gamma)],
	            [0, -sin(gamma), cos(gamma)]])

	Cw = array([[ cos(psi), 0,  -sin(psi)],
	            [      0,   1,          0],
	            [ sin(psi), 0,  cos(psi)]])

	Cv = array([[ cos(theta), sin(theta), 0],
	            [-sin(theta), cos(theta), 0],
                [      0,              0, 1]])

	A = Cy.dot(Cv).dot(Cw)
	return A.transpose()

a = 6378137.0
b = 6356752.3142
h = 170
lambda0 = 37.499780
fi0     = 55.81158

def G(a, b, h, fi0):
	u = 0.7292115e-4
	ge = 9.78049
	fi0 = 55.811527

	q = u ** 2 * a / ge
	e = ((a ** 2 - b ** 2) / a ** 2) ** 0.5
	alpha = (a - b) / a

	betta  = 5/2 * q - alpha - 17/14 * q * alpha
	betta1 = alpha  ** 2 / 8 - 5 / 8 * q * alpha

	g0 = ge * (1 + betta * sin(fi0) ** 2 + betta1 * sin(fi0) ** 2)

	gx = 0
	gy = 0.5 * ge * h / a * (4 * q - e ** 2) * sin(2 * fi0)
	gz = g0 + h / a * (3 * h / a - 2*q * cos(fi0) **2 + e ** 2 * ( 3* sin(fi0) ** 2 - 1) - q * (1 + 6 * sin(fi0) ** 2))

	return array([[gx], [gz], [gy]]) #[gx, gy, gz]


def Read():
	read = []

	with open('data.txt', 'r') as f:
		j = 0

		for line in f:
			LGx, LGy, LGz, AKx, AKy, AKz, pitch, roll, Tyaw, Vn, Ve = [float(i) for i in line.split()]
			read.insert(j, [LGx, LGy, LGz, AKx, AKy, AKz, pitch, roll, Tyaw, Vn, Ve])
			j += 1

	read = array(read)
	read = read.transpose()
	return read


def sys(u, n, V, Lamda, Fi, hh, Lamda_dot, Fi_dot, H_dot, f):#, H_dot):
	dt = 1
	fi = Fi[-1]
	h =170
	g = G(a, b, hh[-1], fi)
	fidot = Fi_dot[-1]
	adot = Lamda_dot[-1]
	Vx = V[0][-1]
	Vy = V[1][-1]
	Vz = V[2][-1]

	Vx_dot = -(2*u*sin(fi) + adot*sin(fi))*Vz + g[0][0] + n[0][0] - fidot*Vy +0.1624588541522454
	Vy_dot = -(2*u*cos(fi) + adot*cos(fi))*Vz - g[1][0] + n[1][0] - fidot*Vx -  0.00859362
	Vz_dot = -(2*u*cos(fi) + adot*sin(fi))*Vy + g[2][0] + n[2][0] + (2*u*sin(fi) + adot*sin(fi))*Vx+0.10509408

	V[0].append(V[0][-1] + Vx_dot * dt)
	V[1].append(V[1][-1] + Vy_dot * dt)
	V[2].append(V[2][-1] + Vz_dot * dt)
	# f.write('Vx= {:.6}    Vy= {:.6}    Vz= {:.6}\n'.format(str(round(V[0][-1],6)), str(round(V[1][-1],6)), str(round(V[2][-1],6))))
	# f.write('Vx\'= {:.6}    Vy\'= {:.6}    Vz\'= {:.6}\n'.format(str(round(Vx_dot,6)), str(round(Vy_dot,6)), str(round(Vz_dot,6))))

	fi_dot =  V[0][-1] / (R + H)
	lamda_dot = V[2][-1] / (R + H) / cos(fi)
	h_dot = V[1][-1]
	Fi_dot.append(fi_dot)
	Lamda_dot.append(lamda_dot)
	H_dot.append(h_dot)

	Fi.append(Fi[-1] + fi_dot * dt)
	Lamda.append(Lamda[-1] + lamda_dot * dt)
	hh.append(hh[-1] + h_dot * dt)
	return Fi[-1], Lamda[-1], hh[-1]


H = 170
Lamda = [RAD(lambda0)]
Fi    = [RAD(fi0)]
HH    = [H]
Lamda_dot = [0]
Fi_dot    = [0]
H_dot     = [0]
V = [[0], [0], [0]]
E = [[], []]

u = 0.7292115e-4
R = 6371000
read = Read()
count = len(read[0])

with open('rez.txt', 'w') as f:
	for i in range(count):
		psi   =read[8][i]
		gamma =read[7][i]
		theta =read[6][i]
		nc = array([[read[3][i]], [read[4][i]], [read[5][i]]])

		A = matr(RAD(psi), RAD(theta), RAD(gamma))
		nH = A.dot(nc)

		fi_temp, lamda_temp, h_temp =sys(u, nH, V, Lamda, Fi, HH, Lamda_dot, Fi_dot, H_dot, f)
		E[0].append(fi_temp - fi0)
		E[1].append(lamda_temp - lambda0)


T = arange(0, count)
fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3)
ax1.plot(T, [fi0] * (count), label='φ0')
ax2.plot(T, [lambda0] * (count), label='λ0')
ax3.plot(T, [H] * (count), label='H0')

Fi = [GRAD(i) for i in Fi]
Lamda = [GRAD(i) for i in Lamda]
print('dFi = ', Fi[-1] - Fi[0])
print('dLamda = ', Lamda[-1] - Lamda[0])
T = arange(0, count+1)
ax1.plot(T, Fi,'-.', label='φ вычисленная')
ax2.plot(T, Lamda,'-.', label='λ вычисленная')
ax3.plot(T, HH,'-.', label='H вычисленная')
ax1.set_ylabel('φ[°]')
ax1.set_xlabel('t, [c]')
ax2.set_ylabel('λ[°]')
ax2.set_xlabel('t, [c]')
ax3.set_ylabel('H[m]')
ax3.set_xlabel('t, [c]')
ax1.grid()
ax2.grid()
ax3.grid()
ax1.legend()
ax2.legend()
ax3.legend()
fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3)
ax1.plot(T, V[0], label='Vx')
ax2.plot(T, V[1], label='Vy')
ax3.plot(T, V[2], label='Vz')
ax1.set_ylabel('Vx[m/c^2]')
ax1.set_xlabel('t, [c]')
ax2.set_ylabel('Vy[m/c^2]')
ax2.set_xlabel('t, [c]')
ax3.set_ylabel('Vz[m/c^2]')
ax3.set_xlabel('t, [c]')
ax1.grid()
ax2.grid()
ax3.grid()
ax1.legend()
ax2.legend()
ax3.legend()

pyplot.show()
