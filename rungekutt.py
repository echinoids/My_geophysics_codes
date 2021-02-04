# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 20:09:56 2020

@author: steph
"""

import numpy as np
import matplotlib.pyplot as plt

#Datos:
R0 = 1.e5
theta = 3.0
Ms = 1.989e30
G = 6.67e-11
sigma_B = 5.67e-8 #Ctte de boltzman

M_p=np.array([4.87,5.97,0.642])*1e24 #[kg]
R_p=np.array([6.05,6.37,3.40])*1e6   #[km]
rho_p=np.array([5.25,5.51,3.90])*1e3 #[kg/m^3]
a_p=np.array([1.08,1.50,2.28])*1e11  #[km]
da_p=np.array([4.49,6.32,7.8])*1e10  #[km]
Cp =1e3                              #[J]

T_ac=100e6                           #Tiempo de acreción aprox [años]
dt=1e3                               #Intervalo de tiempo o paso                     

#Defino condiciones iniciales
M0_p=rho_p*R0**3*4.0*np.pi/3         #[kg]
T0= 300                              #[K]

#Genero las funciones

dM=lambda M,R,Mp,A: A*(Mp-M)*R**2
dT=lambda T,R,M,dM,T0,B,C : B*dM/R + (T0-T)*dM/M + C*(T**4-T0**4)/R

t=np.arange(0,T_ac,dt)*1e-6

R=np.zeros(len(t))
M=np.zeros(len(t))
T=np.zeros(len(t))
    
R[0]=R0
M[0]=rho_p[ip]*4*np.pi*R[0]**3/3
T[0]=T0
    
a=a_p[ip]
da=da_p[ip]
Rp=R_p[ip]
Mp=M_p[ip]
    
vk= np.sqrt(G*Ms/a)
CR = (1+2*theta)*vk*Rp**3/6/a**2/da
    
A=(1+2*theta)*vk/2/a**2/da;print(B0)
B=G*(1+1/2/theta)/Cp
C=-3*sigma_B/Cp/rho_p[ip]