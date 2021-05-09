import numpy as np
np.set_printoptions(precision = 4, linewidth = 150)
import matplotlib.pyplot as plt

f = lambda x,T,L,A,E: -T*(x**3 - L**3)/(6*A*E)

Largo = 1.50
T = -200000
elementos = 3
A = 0.001
E = 210e9

Desplazamientos = np.zeros(elementos + 1)
Fuerzas = np.zeros(elementos + 1)

VectLong = np.linspace(0, Largo, elementos+1)

KGlob = np.zeros((elementos + 1,elementos + 1))
Kloc = np.eye(2)
Kloc[0,1] = -1
Kloc[1,0] = -1

Ftot = T*(VectLong[1]**2 - VectLong[0]**2)/2
Fvec = np.zeros(2)

for i in range(elementos):

    #MDF-COMMENT impecable
	Fuerzas[i] += Ftot/3 + (VectLong[i+1]-VectLong[i])*T*VectLong[i]/2
	Fuerzas[i+1] += Ftot*2/3 + (VectLong[i+1]-VectLong[i])*T*VectLong[i]/2

	K = E*A/(VectLong[i+1]-VectLong[i])

	KGlob[i:i+2,i:i+2] += Kloc*K

print(np.arange(elementos))

#MDF-COMMENT en realidad cuando estas adentro de un paréntesis no hace falta la continuacion de linea.
#MDF-COMMENT a demas por legibliidad conviene escribir un poco distinto:
#MDF-COMMENT Desplazamientos[0:elementos] = np.linalg.solve(KGlob[np.ix_(np.arange(elementos),\
#MDF-COMMENT 	np.arange(elementos))], Fuerzas[np.arange(elementos)])
Desplazamientos[0:elementos] = np.linalg.solve(
        KGlob[np.ix_(np.arange(elementos), np.arange(elementos))],
        Fuerzas[np.arange(elementos)] #MDF-COMMENT y acordate que aca te falta la parte delos desplazamientos
        )
Respuesta = np.dot(KGlob[elementos,np.arange(elementos + 1)], Desplazamientos) - Fuerzas[elementos]

plt.plot(np.linspace(0, Largo, 100), f(np.linspace(0, Largo, 100),T,Largo,A,E),'r')
plt.plot(np.linspace(0, Largo, elementos+1), Desplazamientos)
#MDF-COMMENT y te quedaria ponerle una leyenda para que se entinenda.

#MDF-COMMENT preguntas:
#MDF-COMMENT 1. Como evoluciona cuando cambias el numero de elementos ?
#MDF-COMMENT 2. te falta calcular las tensiones no ?
plt.show()
