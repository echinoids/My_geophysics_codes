# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:08:56 2021

@author: steph
"""
import numpy as np

dat=np.loadtxt('planetas.txt')

#Datos almacenados en vectores

Mp=6.5e24                   #Masa final del planeta [kg]
rhop=5476                   #Densidad del planeta [kg/m^3]
ap=9e10                     #Radio orbital [m]
Rp=0.999*6.37e6             #Radio de cada planeta [m] (limite superior)
R0=0.01*Rp                  #Radio inicial de cada planeta 0.01*Rp (limite inferior)[m]
dap=4.5e8

def croot(x):
    x = abs(x)
    return x**(1./3.)

def acrecion(rhop,ap,dap,Rp,Mp,R0):

    #Constantes
    theta = 3.0              #Parametro de Safronov
    P=(1+2*theta)
    Ms = 1.989e30            #Masa del sol [kg/m^3]
    G = 6.671e-11            #Constante de gravitación universal [Nm^2/kg^2]
    
    vk= np.sqrt(G*Ms/ap)     #Velocidad kepleriana [m/s]
    A=4*np.pi*rhop
    B=(P*vk)/(2*dap*ap**2)
    C=(A/(B*Mp))*croot(3*Mp/A)
    D=croot((3*Mp)/A)
    X=Rp/D
    X0=R0/D
    E=np.sqrt(1/3)
    
    #Tiempo de acreción
    
    #En R(t=t')=Rp
    t1=(-1/3)*np.log(np.abs(1-X))
    t2=(1/6)*np.log(np.abs(X**2+X+1))
    t3=E*np.arctan(2*E*X+E)
    
    #En R(t=0)=0.01Rp
    t10=(-1/3)*np.log(np.abs(1-X0))
    t20=(1/6)*np.log(np.abs(X0**2+X0+1))
    t30=E*np.arctan(2*E*X0+E)
    
    t=C*((t1+t2+t3)-(t10+t20+t30))                  #en segundos
    
    return t/(60*60*24*365.25)                      #en años

#Calculo

t_ac=acrecion(rhop,ap,dap,Rp,Mp,R0)
    
print("El tiempo de acrecion para cada planeta es:",t_ac/1e6,"[Ma]")

#Venus, Tierra, Marte#