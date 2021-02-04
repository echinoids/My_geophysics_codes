# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 02:44:08 2021

@author: steph
"""

import numpy as np
import matplotlib.pyplot as plt
from ConveccionTierraCapaLimite import ConveccionTierraCapaLimite

"""
     Evolucion temporal de la temperatura del manto y del nucleo bajo la
     hipotesis de capa limite
    
     Inputs
    
     R = radio de la Tierra [m]
     Rc = radio del nucleo [m]
     T0 = Temperatura exterior [K]
     Rac = Numero de Rayleigh critico
     C = parametros nucleo (C_c) o manto (C_m) {array de 5x1}
           - C[0] = difusividad termica [m**2/s]
           - C[1] = conductividad termica [W/K/m]
           - C[2] = Coef expansion termica [1/K]
           - C[3] = viscosidad cinematica [m**2/s]
           - C[4] = gravedad [m/s**2]
     n0 = numero de isotopos a t = t0 {array 4x1}
           - n0[0] = isotopos de 238U
           - n0[1] = isotopos de 235U
           - n0[2] = isotopos de 232Th
           - n0[3] = isotopos de 40K
     fm = fraccion de isotopos en el manto
     fc = fraccion de isotopos en el nucleo
     t0 = tiempo cero [a~nos]
     tf = tiempo final [a~nos]
     dt = intervalo de integracion [a~nos]
     Tc0 = Temperatura nucleo al tiempo t0 [K]
     Tm0 = Temperatura manto al tiempo t0 [K]
    
     Outputs
    
     t = tiempo {array (tf/dt + 1)x1}[a~nos]  
     Tc = temperatura nucleo {array (tf/dt + 1)x1}[K]
     Tm = temperatura manto {array (tf/dt + 1)x1}[K]
    """
R=6.371*1e6
Rc=3.48*1e6
T0=270
Rac=1e3
Cc=np.array([6*1e-6,36,1e-5,1,7])
Cm=np.array([1e-6,6,1.5*1e-5,1e18,10])
fm=0.7
fc=0.
t0=0
tf=4*1e9
dt=1e7
#Tc0=4500 #Nucleo original
#Tc0=5000 #Nucleo caliente
Tc0=3000 #Nucleo frio
Tm0=3000

#Calculo del numero de isotopos para t=0

dat=np.loadtxt('tabla.txt')

#Constantes importantes

Na=6.02*1e23                       #Número de avogadro

#----------------------------------------------------------------------
#Caso para un planeta sin estratificación
#----------------------------------------------------------------------

#Número de partículas - Energía de cada isotopo y total - Fuente de calor

#Defino vectores y unidades [m.k.s]

masa=dat[:,2]     #Masa total [kg]
As=dat[:,3]       #Masa atómica de los elementos [kg]
masa_m=masa*0.7   #Masa del manto [kg]
lamb=dat[:,1]  #Constante de decaimiento[1/s]
t=(4*1e9)*365*24*3600 #Tiempo en el pasado [s]

N=np.zeros(len(masa))    #N° de particulas de la tierra actualmente
n0=np.zeros(len(masa))   #N° de particulas de la tierra hace 4x10^9 años

for i in np.arange(len(masa)):
    N[i]=(masa[i]/As[i])*Na
    n0[i]=N[i]*np.exp(lamb[i]*t)

[t,Tc,Tm] =  ConveccionTierraCapaLimite(R,Rc,T0,Rac,Cc,Cm,n0,fm,fc,t0,tf,dt,Tc0,Tm0)

datos=[Tc,Tm,Tm/Tc]
name=['Tc','Tm','Tm/Tc']
name2=['Tc [K]','Tm [K]','Tm/Tc [K]']
axy0=[np.min(Tm),np.max(Tc)]
axy1=[np.min(Tm),np.max(Tc)]
axy2=[0,1]
axx0=[0,tf]
axx1=[0,tf]
axx2=[0,tf]
axisy=[axy0,axy1,axy2]
axisx=[axx0,axx1,axx2]

fig, axes = plt.subplots(1, len(datos),figsize=(10,5))
fig.subplots_adjust(hspace=1, wspace=0.4)

for ax,axisx,axisy,datos,name,name2 in zip(axes.flat,axisx,axisy
                                           ,datos,name,name2):
    ax.plot(t,datos,'k')
    ax.set_title(name)
    ax.set_xlabel("t [anos]")
    ax.set_ylabel(name2)
    ax.ticklabel_format(useOffset=False)
    ax.set_xlim(axisx)
    ax.set_ylim(axisy)
plt.show()