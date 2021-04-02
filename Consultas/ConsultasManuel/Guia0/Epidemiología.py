import numpy as np
import matplotlib.pyplot as plt

S0 = 0.9
I0 = 0.1
R0 = 0.0

def funcion (i):
       #MDF-COMMENT lo que no me gusta de este approach es que necesitas 3 llamadas a la función ,
       #MDF-COMMENT cuando podrías hacerlo con una sola
       #MDF-COMMENT y por otro lado, los valores del paso anterior (S,I,R) tendrían que ser los verdaderos 
       #MDF-COMMENT argumentos de tu función.
	dS = lambda t,S,I,R: -0.35*I*S
	dI = lambda t,S,I,R: 0.35*I*S - 0.1*I
	dR = lambda t,S,I,R: 0.1*I
	Der = [dS, dI, dR]

	return Der[i]

V#MDF-COMMENT a ver si me sale.
def funcionv(Y):
    """
    Parámetros
    ==========
    Y = [S, I, R]
    """
    dY = [ 
            -0.35*Y[1]*Y[0],
            0.35*Y[1]*Y[0]-0.1*Y[1],
            0.1*Y[1]
        ]
    return dY



t = np.linspace(0, 100, 1000)
k = np.zeros(4)
Y = np.zeros((3,len(t)))
Y[0,0] = S0
Y[1,0] = I0
Y[2,0] = R0

Espaciado = t[1] - t[0]

for i in range(1,len(t)):
	for j in range(0,3):
		G = funcion(j)
		k[0] = G(t[i-1], Y[0,i-1], Y[1,i-1], Y[2,i-1])
		k[1] = G(t[i-1] + Espaciado/2 , Y[0,i-1] + k[0]*Espaciado/2, Y[1,i-1] + k[0]*Espaciado/2, Y[2,i-1] + k[0]*Espaciado/2)
		k[2] = G(t[i-1] + Espaciado/2 , Y[0,i-1] + k[1]*Espaciado/2, Y[1,i-1] + k[1]*Espaciado/2, Y[2,i-1] + k[1]*Espaciado/2)
		k[3] = G(t[i-1] + Espaciado, Y[0,i-1] +k[2]*Espaciado, Y[1,i-1] +k[2]*Espaciado, Y[2,i-1] +k[2]*Espaciado)
		Y[j,i] = Y[j,i-1] + 1/6*(k[0] + 2*k[1] + 2*k[2] + k[3])*Espaciado

plt.plot(t,Y[0,:])
plt.plot(t,Y[1,:])
plt.plot(t,Y[2,:])
plt.show()
