import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt
import pdb #MDF-COMMENT for debugging
#MDF-COMMENT para que me funcione con la estructura de archivos que pusiste 
#MDF-COMMENT hay que ponerle el directorio actual.
import os

F, u, vinculos = np.loadtxt(
        os.path.join('Datos Ejercicio2', 'DatosFuV.dat'),
        delimiter = ',\t',
        unpack = True
        )
        #MDF-COMMENT"DatosFuV.dat", delimiter=',	', unpack = True)
E, A = np.loadtxt(
        os.path.join('Datos Ejercicio2',"DatosEA.dat"),
        delimiter=',\t',
        unpack = True)
# MatrizNodos = np.loadtxt("MatrizNodos.dat", delimiter=',	')
MatrizNodos = np.loadtxt(
       os.path.join('Datos Ejercicio2', "MatrizNodos.dat"),
       delimiter=',\t',
       dtype=int)
MatrizConectividad = np.loadtxt(
       os.path.join('Datos Ejercicio2', "MatrizConectividad.dat"),
       delimiter=',\t',
       dtype=int)
#MDF-COMMENT dimensiones = int(np.shape(MatrizNodos)[1])
#MDF-COMMENT nodos = int(np.shape(MatrizNodos)[0])
nodos, dimensiones = MatrizNodos.shape
#MDF-COMMENT elementos = int(np.shape(MatrizConectividad)[0])
elementos = MatrizConectividad.shape[0]
Tensiones = np.zeros(elementos)

#MDF-COMMENT MatrizConectividad = MatrizConectividad.astype(int)
#MDF-COMMENT no hace falta porque ya pusimos para leer int directamente

Klocal = np.zeros((dimensiones, dimensiones))
KGlobal = np.zeros((dimensiones*nodos, dimensiones*nodos))
k = []
#LargoV = np.zeros((elementos,2))
angulo = []

for i in range(elementos):
#MDF-COMMENT mind blown:
#MDF-COMMENT for e, MCL in enumerate(MatrizConectividad):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
        #MDF-COMMENT ... = MatrizNodos[MCL[0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
        #MDF-COMMENT ... = MatrizNodos[MCL[1],:]
	N1, N2 = MatrizConectividad[i,:]
	Largo = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	k.append(E[i]*A[i]/Largo)
	#k = E[i]*A[i]/Largo
	angulo.append(np.arctan2(y2-y1, x2-x1))
	sen = np.sin(angulo[i])
	cos = np.cos(angulo[i])

	#LargoV[i,0] = Largo

	Klocal[0,0] = cos**2
	Klocal[1,1] = sen**2
	Klocal[0,1] = cos*sen
	Klocal[1,0] = cos*sen
	Klocal = Klocal*k[i]

	#print(Klocal)
        #MDF-COMMENT empecable
        #MDF-COMMENT sin embargo, te recomiendo que trates de implementar el esquema 
        #MDF-COMMENT general que hablamos, porque no vas a poder hacer esto cuando tengamos matrices de más dimensiones.
        #MDF-COMMENT es cierto que lo del rango funciona perfecto, me tengo que acordar para el año que viene,
        #MDF-COMMENT aunque a mí personalmente me marea un poco que hay que ponerle uno más que el rango real
	KGlobal[2*N1:(2*N1 + 2), 2*N1:(2*N1 + 2)] += Klocal
	KGlobal[2*N2:(2*N2 + 2), 2*N2:(2*N2 + 2)] += Klocal
	KGlobal[2*N1:(2*N1 + 2), 2*N2:(2*N2 + 2)] += -Klocal
	KGlobal[2*N2:(2*N2 + 2), 2*N1:(2*N1 + 2)] += -Klocal

	plt.plot([x1, x2], [y1, y2], 'k', '-')
	plt.text((x2+x1)/2,(y2+y1)/2, i, fontsize = 10, color = 'red')


SubF = []
R = []
S = []

for i in range(dimensiones*nodos):
	if vinculos[i] == 1:
		SubF.append(F[i])
		R.append(i)
	else:
		S.append(i)

#MDF-COMMENT en este caso está bien, pero acordate que en el caso general si tenes vínculos de desplazamietnto
#MDF-COMMENT distintos de cero vas a tener que agregar el término en la línea esta
u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R])
#MDF-COMMENT u[R] = np.linalg.solve(KGlobal[np.ix_(R,R)], F[R] - K[np.ix_(R, S)].dot(u[s]))
#MDF-COMMENT en este otro caso tampoco necesitas el np.ix_:
#MDF-COMMENT F[S] = KGlobal[S,:].dot(u)
F[S] = np.dot(KGlobal[np.ix_(S,np.arange(dimensiones*nodos))],u)

U = []
Floc = np.zeros(4)
Kext = np.eye(4)

#MDF-COMMENT 
for i in range(elementos):
	x1, y1 = MatrizNodos[MatrizConectividad[i,0],:]
	x2, y2 = MatrizNodos[MatrizConectividad[i,1],:]
	N1, N2 = MatrizConectividad[i,:]

	x1 += u[2*N1]
	x2 += u[2*N2]
	y1 += u[2*N1+1]
	y2 += u[2*N2+1]

	U = [u[2*N1], u[2*N1+1], u[2*N2], u[2*N2+1]]

	#Largo = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	#LargoV[i,1] = Largo

	#Tensiones[i] = (np.dot([F[2*N2] - F[2*N1], F[2*N2 + 1] - F[2*N1 + 1]], [(x2-x1),(y2-y1)]))/A[i]
	x1 += u[2*N1]*149
	x2 += u[2*N2]*149
	y1 += u[2*N1+1]*149
	y2 += u[2*N2+1]*149

	sen = np.sin(angulo[i])
	cos = np.cos(angulo[i])

	L = np.arange(4)

        #MDF-COMMENT aca no entiendo
	Kext[L,L] = cos
	Kext[0,1] = sen
	Kext[1,0] = -sen
	Kext[2,3] = sen
	Kext[3,2] = -sen

	Floc = np.dot(k[i]*Kext, U)
	Tensiones[i] = (Floc[2] - Floc[0])/A[i]


	
	plt.plot([x1, x2], [y1, y2], 'b', '-')



for i in range(nodos):
	plt.text(MatrizNodos[i,0],MatrizNodos[i,1], i, fontsize = 14)

print(Tensiones)
#print(F)
#print(u)
#print(LargoV)
plt.legend(["Antes", "Después"])
#MDF-COMMENT impecable

#MDF-COMMENT y ya que estas guardá la figura, no ?
plt.savefig('Puente_Manuel.pdf')
# plt.show()

plt.close()
# pero mas facil:
NodosDESP = MatrizNodos + u.reshape((nodos,2), order='C')*np.array([500, 100])
#MDF-COMMENT                arma [[ux0, ux1], [ux1, ux2], etc ]  y multiplica cada columna por una escala distinta
plt.figure()
for e , elem in enumerate(MatrizConectividad):
    before=plt.plot(MatrizNodos[elem, 0], MatrizNodos[elem, 1], 'ko-', ms = 10 , label='sin desplazar')[0]
    after =plt.plot(NodosDESP[elem, 0], NodosDESP[elem, 1], 'ro-', ms = 10 ,  label='desplazado' )[0]
plt.legend([before, after], [before.get_label(), after.get_label()] )
plt.show()
