# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 17:11:02 2021

@author: steph
"""
import numpy as np

def croot(x):
    x = abs(x)
    return x**(1./3.)

dat=np.loadtxt('planetas.txt')

#Datos almacenados en vectores

Mp=dat[:,0]              #Masa final de cada planeta [kg]
Rp=dat[:,1]              #Radio de cada planeta [m]
rhop=dat[:,2]            #Densidad del planeta=densidad de los planetesimales [kg/m^3]

#Constantes

theta = 3.0              #Parametro de Safronov
G = 6.671e-11            #Constante de gravitación universal [Nm^2/kg^2]
Cp =1e3                  #Capacidad específica [J]
gamma=0.1

dT=np.zeros(len(dat[:,0]))

for i in range(len(dat[:,0])):
    
    A=1/(Cp*(1+gamma))
    B=(G*gamma**(2/3))/6
    C=(1+gamma**(5/3)-(1+gamma)**(5/3))*(3*G/5)
    
    dT[i]=A*(B*(Mp[i]/Rp[i])-C*(Mp[i]**(2/3))*croot(4*np.pi*rhop[i]/3))


print("La diferencia de temperatura es:",dT, "[K]")