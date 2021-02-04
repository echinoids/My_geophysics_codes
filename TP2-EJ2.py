# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 21:01:42 2020

@author: steph
"""

#Ejercicio 2

import numpy as np
import matplotlib.pyplot as plt

dat=np.loadtxt('planetas.txt')

#Datos almacenados en vectores

Mp=dat[:,0]                 #Masa de cada planeta [kg]
Rp=dat[:,1]                 #Radio de cada planeta [m]
rhop=dat[:,2]               #Densidad del planeta [kg/m^3]
ap=dat[:,3]                 #Radio del anillo de cada planeta [m]
dap=dat[:,4]                #Ancho del anillo de cada planeta [m]

#Condiciones iniciales
R0 = 1e5                    #Radio [m]
T0=300                      #Temperatura inicial [K]
T_ac=1200e5                 #Tomamos un valor de referencia ya que los del EJ1 son muy chicos.
dt=1e3                      #Paso de integración [años]
Nit=np.int(T_ac/dt)         #Numero de iteraciones
t=np.arange(0,T_ac,dt)*1e-6 #Vector de tiempo [Maños]

#Genero las funciones (Es una forma alternativa a la usual def function ():)

#dR=lambda R,Rp,ctte1 : ctte1*(1.0-(R/Rp)**3)
#dM=lambda M,R,Mp,ctte2 : ctte2*(Mp-M)*R**2
#dT=lambda T,R,M,dM,T0,ctte3,ctte4 : ctte3*dM/R + (T0-T)*dM/M - (ctte4/R)*(T**4-T0**4)

def f(R,Rp,M,Mp,T,T0,ap,dap,rhop):
    
    theta = 3.0              #Parametro de Safronov
    Ms = 1.989e30            #Masa del sol [kg/m^3]
    G = 6.671e-11            #Constante de gravitación universal [Nm^2/kg^2]
    year = 365*24*60*60      #Ctte para llevar de segundos a años [s/año]
    sigma = 5.67e-8          #Constante de Boltzman [W/m^2K^4]
    Cp =1e3                  #Capacidad específica [J]
    
    #Defino las constantes
    vk= np.sqrt(G*Ms/ap)*year
    ctte1 = (1+2*theta)*vk*Rp**3/(6*ap**2*dap)
    ctte2 = ((1+2*theta)*vk)/(2*ap**2*dap)
    ctte3 = (G/Cp)*(1+1/(2*theta))
    ctte4 = (3*sigma*year)/(rhop*Cp)
    
    #Defino incrementos
    dR=ctte1*(1-(R/Rp)**3)
    dM=ctte2*(Mp-M)*R**2
    dT=ctte3*dM/R + (T0-T)*dM/M - ctte4*(T**4-T0**4)/R

    return dR,dM,dT

#Defino matrices donde voy a almacenar la información

R=np.zeros((Nit,len(dat[:,0])))
M=np.zeros((Nit,len(dat[:,0])))
T=np.zeros((Nit,len(dat[:,0])))

#Resuelvo la EDO por medio del metodo de Euler

for i in range(len(dat[:,0])):
    
    #Defino que la primera muestra de los vectores anteriores sean
    #las condiciones iniciales
    R[0,i]=R0
    M[0,i]=rhop[i]*(4/3)*np.pi*R0**3
    T[0,i]=T0
    
    for j in range(Nit-1):    #Debe ser Nit-1 porque hemos colocado el valor en la posición [0]
        
        dR,dM,dT=f(R[j,i],Rp[i],M[j,i],Mp[i],T[j,i],T0,ap[i],dap[i],rhop[i])
        
        R[j+1,i]=R[j,i]+dt*dR
        M[j+1,i]=M[j,i]+dt*dM
        T[j+1,i]=T[j,i]+dt*dT

#Graficos de los resultados

name=['Venus','Tierra','Marte']
ylabel=['Radio [Km]','Masa [Kg]','Temperatura [K]']

Av=R[:,0];Bv=M[:,0];Cv=T[:,0];At=R[:,1];Bt=M[:,1];Ct=T[:,1];Am=R[:,2];Bm=M[:,2];Cm=T[:,2]

V=[Av*1e-3,Bv,Cv]
T=[At*1e-3,Bt,Ct]
Ma=[Am*1e-3,Bm,Cm]

fig, axes = plt.subplots(len(dat[:,0]),1,figsize=(8,9))
fig.subplots_adjust(hspace=1, wspace=0.2)

for ax,V,T,Ma,ylabel in zip(axes.flat,V,T,Ma,ylabel):
    ax.plot(t,V,label='Venus')
    ax.plot(t,T,label='Tierra')
    ax.plot(t,Ma,label='Marte')
    ax.set_xlabel('t [Ma]')
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    
plt.show()

