from math import pi, cos, sin, sqrt, atan
from numpy import array, dot, transpose, empty


def RAD(grad):
	return grad / 180 * pi


def A(A0, U, t):
	return RAD(A0 + U * t)


def G(Re, H, fi, lamb):
	zg = (Re + H) * sin(RAD(fi))
	yg = (Re + H) * cos(RAD(fi)) * sin(RAD(lamb))
	xg = (Re + H) * cos(RAD(fi)) * cos(RAD(lamb))
	return array([[xg], [yg], [zg]])


def modylVp(Vn):
	return sqrt(Vn[0][0]**2 + Vn[2][0]**2)


def Vg(Vxg, Vyg, Vzg):
	return array([[Vxg], [Vyg], [Vzg]])


def all(i):
	fi = (90 - 180 / 22 * i)*pi/180
	lambd = (360 - 360 / 22 * i)*pi/180
	H = 3000 / 22 * i
	A0 = 180 + 180 / 22 * i
	U = 15
	t = 24 - 24 / 22 * i
	Vxg = -4 + 8/22 * i
	Vyg = -1 + 2/22 * i
	Vzg = -3 + 6/22 * i
	print("Vxg =", Vxg)
	print("Vyg =", Vyg)
	print("Vzg =", Vzg)

	Ag1 = array([[cos(lambd), sin(lambd),  0],
              [-sin(lambd), cos(lambd),	0],
              [0,          0,  1]])

	A12 = array([[cos(fi), 0, -sin(fi)],
              [0, 1,        0],
              [sin(fi), 0, cos(fi)]])

	Ap = array([[0, 0, 1],
             [1, 0, 0],
             [0, 1, 0]])

	matrica1 = Ap.dot(A12)
	matrica2 = matrica1.dot(Ag1)
	'''matrica2_t = matrica2.transpose()'''
	print('1:\n', matrica1)
	print('2:\n', matrica2)

	Vg_1 = Vg(Vxg, Vyg, Vzg)

	Vn = matrica2.dot(Vg_1)
	print('Vn =\n', Vn)
	print('\nVn{1,1}')
	print(Vn[0][0])
	print('\nVn{1,2}')
	print(Vn[1][0])
	print('\nVn{1,3}')
	print(Vn[2][0])
	print('1: |Vn| = ', modylVp(Vn))
	#-------------------
	if Vn[1][0] > 0:
		print('2: Up')
	elif Vn[1][0] < 0:
		print('2: Down')
	else:
		print('2: Aflat')

	#_________________________
	fi_p = atan(Vn[2][0] / Vn[0][0])
	print('3: fi_p = ', fi_p * 180 / pi)


all(8)
