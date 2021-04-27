import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt

F, u, vinculos = np.loadtxt("DatosFuV.dat", delimiter=',	', unpack = True)
E, A = np.loadtxt("DatosEA.dat", delimiter=',	', unpack = True)
MatrizNodos = np.loadtxt("MatrizNodos.dat", delimiter=',	')
MatrizConectividad = np.loadtxt("MatrizConectividad.dat", delimiter=',	')

dimensiones = int(np.shape(MatrizNodos)[1])
nodos = int(np.shape(MatrizNodos)[0])
elementos = int(np.shape(MatrizConectividad)[0])

MatrizConectividad = MatrizConectividad.astype(int)

Klocal = np.zeros((dimensiones, dimensiones))
KGlobal = np.zeros((dimensiones*nodos, dimensiones*nodos))

for i in range(elementos):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
	N1, N2 = MatrizConectividad[i,:]
	Largo = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	k = E[i]*A[i]/Largo
	angulo = np.arctan2(y2-y1, x2-x1)
	sen = np.sin(angulo)
	cos = np.cos(angulo)

	Klocal[0,0] = cos**2
	Klocal[1,1] = sen**2
	Klocal[0,1] = cos*sen
	Klocal[1,0] = cos*sen
	Klocal = Klocal*k

	KGlobal[2*N1:(2*N1 + 2), 2*N1:(2*N1 + 2)] += Klocal
	KGlobal[2*N2:(2*N2 + 2), 2*N2:(2*N2 + 2)] += Klocal
	KGlobal[2*N1:(2*N1 + 2), 2*N2:(2*N2 + 2)] += -Klocal
	KGlobal[2*N2:(2*N2 + 2), 2*N1:(2*N1 + 2)] += -Klocal

	plt.plot(np.linspace(x1, x2, 2), np.linspace(y1, y2, 2), 'k', '-')

SubF = []
I = []

for i in range(dimensiones*nodos):
	if vinculos[i] == 1:
		SubF.append(F[i])
		I.append(i)

I = np.array(I)
SubU = np.zeros(len(I))
SubF = np.ravel(SubF)
SubK = KGlobal[np.ix_(I,I)]

SubU = np.linalg.solve(SubK, SubF)

for i in range(len(I)):
	u[I[i]] = SubU[i]

F = np.dot(KGlobal,u)

for i in range(elementos):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
	x1 += u[2*MatrizConectividad[i,0]]
	x2 += u[2*MatrizConectividad[i,1]]
	y1 += u[2*MatrizConectividad[i,0]+1]
	y2 += u[2*MatrizConectividad[i,1]+1]

	plt.plot(np.linspace(x1, x2, 2), np.linspace(y1, y2, 2), 'b', '--')

for i in range(nodos):
	plt.text(MatrizNodos[i,0],MatrizNodos[i,1], i, fontsize = 14)

plt.legend(["Antes", "Después"])
plt.show()