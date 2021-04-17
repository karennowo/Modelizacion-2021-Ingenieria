import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

Deltat = 1
TiempoF = 30
Longitud = 10
Deltax = 0.1
k = 0.835
Lambda = k*Deltat/(Deltax**2)
xsize = int(np.floor(Longitud/Deltax))
tsize = int(np.floor(TiempoF/Deltat))
T0 = 100
T1 = 50

"""
MallaT = np.zeros((tsize,xsize))
MallaT[:,0] = T0
MallaT[:,-1] = T1
# MallaL = np.eye(xsize)*(1-2*Lambda)

for i in range(1,tsize):
	for j in range(1,xsize-1):
		MallaT[i,j] = MallaT[i-1,j] + Lambda*(MallaT[i-1,j+1] - 2* MallaT[i-1,j] + MallaT[i-1,j-1])

X = np.linspace(0, Longitud, xsize)

fig = plt.figure()

def animacion (i):
	plt.plot(X, MallaT[i,:])

grafica = animation.FuncAnimation(fig, animacion, range(0,tsize))
plt.show()

"""

def animacion2 (i,X,T,A):
	axe[A].plot(X, T[i])

def animacion3 (i,E,A):
	axe[A].semilogy(i,E[i],'k*')

def animacion2n (i,X,T,A):
	axe[A].plot(X, T[i])

def animacion3n (i,E,A):
	axe[A].semilogy(i, E[i],'r*')

X = np.linspace(0, Longitud, xsize)

TempI = np.zeros(xsize)
TempI2 = np.zeros(xsize)
TempIN = np.zeros(xsize)
TempIN2 = np.zeros(xsize)
TempI[0] = T0
TempI[-1] = T1
TempIN[0] = T0
TempIN[-1] = T1
count = 0
T = [TempI]
Tn = [TempIN]
E = [1]
EN = [1]

Mallatiempo = np.eye(xsize)
Mallatiempo2 = np.eye(xsize)
Mallatiempo3 = np.eye(xsize)

for i in range(1,xsize-1):
	Mallatiempo[i,i] =  (1 + 2*Lambda)
	Mallatiempo[i,i-1] = - Lambda
	Mallatiempo[i,i+1] = - Lambda

	Mallatiempo2[i,i] =  2*(1 - Lambda)
	Mallatiempo2[i,i-1] = Lambda
	Mallatiempo2[i,i+1] = Lambda

	Mallatiempo3[i,i] =  2*(1 + Lambda)
	Mallatiempo3[i,i-1] = -Lambda
	Mallatiempo3[i,i+1] = -Lambda


while count < 1000:
	TempI2 = np.linalg.solve(Mallatiempo, TempI)
	E.append(np.dot((TempI - TempI2),(TempI - TempI2)))

	if E[-1] < 1e-5:
		break

	TempI = TempI2
	T.append(TempI)

	count += 1

print(count)


count2 = 0

while count2 < 1000:
	TempIN2 = np.linalg.solve(Mallatiempo3, np.dot(Mallatiempo2,TempIN))
	EN.append(np.dot((TempIN - TempIN2),(TempIN - TempIN2)))

	if EN[-1] < 1e-5:
		break

	TempIN = TempIN2
	Tn.append(TempIN)

	count2 += 1

print(count2)

fig2,axe = plt.subplots(1,2,figsize = (20,10))
axe[1].set_xlim([0,len(EN)])
axe[1].set_ylim([min(EN[1:]),max(EN[1:])])

axe[0].set_ylabel("Temperatura")
axe[0].set_xlabel("Posición")

axe[1].set_ylabel("Error ($|T_{i} - T_{i-1}|$)")
axe[1].set_xlabel("Iteración")

graf3 = animation.FuncAnimation(fig2, animacion2n, range(len(Tn)), fargs = (X,Tn,0), repeat = False, interval = 1)
graf4 = animation.FuncAnimation(fig2, animacion3n, range(len(EN)), fargs = (EN,1), repeat = False, interval = 1)
plt.show()


plt.rc('axes',labelsize = 10)
fig,axe = plt.subplots(1,2,figsize = (10,5))

axe[0].set_ylabel("Error ($|T_{i} - T_{i-1}|$)")
axe[0].set_xlabel("Iteración")
axe[1].set_ylabel("Error C-N ($|T_{i} - T_{i-1}|$)")
axe[1].set_xlabel("Iteración")

axe[0].set_xlim([0,len(E)])
axe[0].set_ylim([min(E[1:]),max(E[1:])])
axe[1].set_xlim([0,len(EN)])
axe[1].set_ylim([min(E[1:]),max(E[1:])])

graf1 = animation.FuncAnimation(fig, animacion3, range(len(E)), fargs = (E,0), repeat = False, interval = 1)
graf2 = animation.FuncAnimation(fig, animacion3n, range(len(EN)), fargs = (EN,1), repeat = False, interval = 1)

plt.show()



"""

plt.rc('axes',labelsize = 20)
fig,axe = plt.subplots(1,2,figsize = (20,10))
axe[1].set_xlim([0,len(E)])
axe[1].set_ylim([min(E[1:]),max(E[1:])])

axe[0].set_ylabel("Temperatura")
axe[0].set_xlabel("Posición")

axe[1].set_ylabel("Error ($|T_{i} - T_{i-1}|$)")
axe[1].set_xlabel("Iteración")

graf1 = animation.FuncAnimation(fig, animacion2, range(len(T)), fargs = (X,T), repeat = False, interval = 10)
graf2 = animation.FuncAnimation(fig, animacion3, range(len(E)), fargs = (E,), repeat = False, interval = 10)
plt.show()

"""


