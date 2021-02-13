# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from numpy import linalg as LA

#==================================================================
#Función de matriz de convolución
#==================================================================
def convmtx (w, N):
    
    Nw = len(w)
    Nd = Nw +N -1
    A  = np.zeros(shape=(Nd,N))
    n=0
    for j in range(0,N):
        for i in range(0,Nw):
            A[i+n,j] = w[i]
        n=n+1
        
    return A

#==================================================================
#Aplicación del método de deconvolución con Filtro Wiener
#==================================================================

#Defino los parámetros de la ondícula y del filtro

#w=np.array([1,1/2])
#w=np.array([1,2])
w=np.array([6,7,2])
#w=np.array([3,-2,1])
#w=np.array([1,-2,3])
#w = np.array([36,-2,18,2])  #Ondiculas
#w = np.array([36,-2,18,-25])

#Defino los tamaños

No=len(w)        #Largo de la ondicula
Nf=4                    #Largo del filtro conformador
Nd= No + Nf -1          #Largo de la salida deseada 
                        #o convolución con el filtro
                        
#Defino la salida deseada

d=np.zeros(Nd)
d[0]=1

# =================================================================
# Inversa con filtro Wiener
# =================================================================
# A f = d => A'A f = A'd ==> f = inv(A'A) (A'd)

#Calculo la matriz de convolución
A= convmtx(w,Nf)

#Calculo la matriz de autocorrelación de la ondícula
R= np.dot(np.transpose(A),A)

# =================================================================
# Sin preblanqueo
# =================================================================

#Hallo la matriz inversa de la matriz de autocorrelación
Rinv=inv(R)

#Hallo la correlación entre la ondícula y la salida deseada
g=np.dot(np.transpose(A),d)

#Hallo los coeficientes del filtro Wiener
f_wien=np.dot(Rinv,g)

# =================================================================
# Con preblanqueo
# =================================================================
I=np.eye(Nf)

#Parámetro de regularización - Trade off - Preblanqueo
#mu=0.01
#mu=0.1
#mu=1
mu=0

Rp=R+mu*I

Rinvp=inv(Rp)

f_wienp=np.dot(Rinvp,g)

# =================================================================
# Verificación del resultado
# =================================================================

#Defino la salida deseada aproximada
d_aprox=np.convolve(w,f_wienp,mode='full')

#Defino el vector de error y el error medio cuadrático
error= d_aprox[0:len(d)]-d
e_rms=LA.norm(error)

# =================================================================
# Autocorrelaciones
# =================================================================

#Esta parte es solo para el ejercicio 2

A=np.array([3,-2,1])
B=np.array([1,-2,3])

Acorr=np.correlate(A,A,mode='full')
Bcorr=np.correlate(B,B,mode='full')

fig = plt.figure(figsize=(20, 10))

ax = fig.add_subplot(2, 2, 1)
ax.stem(A,use_line_collection='True')
ax.set_title('Ondicula original')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 2, 2)
ax.stem(Acorr,use_line_collection='True')
ax.set_title('Autocorrelacion')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 2, 3)
ax.stem(B,use_line_collection='True')
ax.set_title('Ondicula original')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 2, 4)
ax.stem(Bcorr,use_line_collection='True')
ax.set_title('Autocorrelacion')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

plt.show()

# =================================================================
# Energía parcial o acumulada en función del tiempo
# =================================================================
A=np.array([3,-2,1])
B=np.array([1,-2,3])

def ac_energy(w):
    E=np.zeros(len(w))
    e=0
    for i in range(len(w)):
        e=e+w[i]**2
        E[i]=e
    return E

A_energy=ac_energy(A)
B_energy=ac_energy(B)

# =================================================================
# Gráficos de las energías acumuladas
# =================================================================
fig = plt.figure(figsize=(20, 10))

ax = fig.add_subplot(2, 1, 1)
ax.stem(A_energy,use_line_collection='True')
ax.set_title('Energía de A')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 1, 2)
ax.stem(B_energy,use_line_collection='True')
ax.set_title('Energía de B')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

plt.show()

# =================================================================
# Gráficos
# =================================================================

#Gráfico de los datos anteriores
fig = plt.figure(figsize=(20, 10))

ax = fig.add_subplot(2, 3, 1)
ax.stem(w,use_line_collection='True')
ax.set_title('Ondicula original')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 3, 2)
ax.stem(d,use_line_collection='True')
ax.set_title('Salida deseada')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 3, 3)
ax.stem(f_wienp,use_line_collection='True')
ax.set_title('Filtro Wiener')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

ax = fig.add_subplot(2, 3, 4)
ax.stem(d_aprox,use_line_collection='True')
ax.set_title('w*f_wien')
ax.set_xlabel("Amplitud")
ax.set_ylabel("Muestras")

ax = fig.add_subplot(2, 3, 5)
ax.stem(error,use_line_collection='True')
ax.set_title('Vector de errores')
ax.set_ylabel("Amplitud")
ax.set_xlabel("Muestras")

plt.show()