# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 16:45:47 2021

@author: steph
"""

#Importo librerias

import numpy as np
import matplotlib.pyplot as plt

dat=np.loadtxt('tabla.txt')

#Constantes importantes

Rn=3480*1e3                        #Radio del núcleo [m]
R=6370*1e3                         #Radio terrestre [m]
Na=6.02*1e23                       #Número de avogadro
Vn= (4./3.)*np.pi*Rn**3            #Volumen del núcleo [m^3]
Vm= (4./3.)*np.pi*(R**3-Rn**3)     #Volumen del manto [m^3]
Vt= (4./3.)*np.pi*R**3             #Volumen de la tierra [m^3]
To=270 #Temperatura inicial [K]
kn=36 #Constante de resistividad del nucleo [W/mK]
km=6  #Constante de resistividad del manto [W/mK]
k=8   #Constante de resistividad de la tierra [W/mK]
sigma=5.67*1e-8 #Constante de Boltzmann [W/m^2K^4]
eps=1
t0=4*1e9
t=t0*365*24*3600 #Tiempo en el pasado [s]

#----------------------------------------------------------------------
#Caso para un planeta con estratificación
#----------------------------------------------------------------------

#Defino vectores y unidades m.k.s

E=dat[:,0]     #Energía [MeV]
lamb=dat[:,1]  #Constante de decaimiento[1/s]
masa=dat[:,2]  #Masa total [kg]
As=dat[:,3]    #Masa atómica de los elementos [kg]

E=E*1.60218*1e-13 #Energía [J]
masa_n=masa*0.3   #Masa del nucleo [kg]
masa_m=masa*0.7   #Masa del manto [kg]

#Número de partículas - Energía de cada isotopo y total - Fuente de calor

Nn=np.zeros(len(masa_n))   #N° de particulas del nucleo
Nm=np.zeros(len(masa_m))   #N° de particulas del manto
en=np.zeros(len(masa_n))   #Energia de cada isotopo del nucleo
em=np.zeros(len(masa_m))   #Energia de cada isotopo del manto
en_sum=np.zeros(1)         #Energia total del nucleo
em_sum=np.zeros(1)         #Energia total del manto

for i in np.arange(len(masa_n)):
    Nn[i]=(masa_n[i]/As[i])*Na
    Nm[i]=(masa_m[i]/As[i])*Na
    en[i]=E[i]*lamb[i]*Nn[i]
    en_sum=en_sum+en[i]
    em[i]=E[i]*lamb[i]*Nm[i]
    em_sum=em_sum+em[i]
    an=en_sum/Vn                  #Fuente de calor del nucleo[J/m^3s]
    am=em_sum/Vm                  #Fuente de calor del manto [J/m^3s]

#Graficos

#Perfil de calor

plt.figure(figsize=(10,6),facecolor='w')

n=100
r=np.arange(0,R,R/n)                #[m] Posiciones en dirección radial
r[0]=1e-10                          #porque en algun momento divido por r

q=np.zeros(len(r))                  #[J/m²s] Calor por unidad de area y tiempo
for j in np.arange(len(r)):
  if r[j]<Rn:
    q[j]=an*r[j]/3.
    plt.plot(r[j]*1e-3,q[j],'k.')
  else:
    q[j]=am*r[j]/3.+(an-am)/3.*Rn**3/(r[j]**2)
    plt.plot(r[j]*1e-3,q[j],'r.')

plt.xlabel('Radio terrestre [km]')
plt.ylabel('q(r) [ J/m²s ]')

plt.show()

#Perfil de temperatura

Ts=(am*R/3.*sigma*eps+((an-am)*Rn**3)/(3.*sigma*eps*R**2)+To**4)**(1/4)

plt.figure(figsize=(10,6),facecolor='w')

n=100
r=np.arange(0,R,R/n)                #[m] Posiciones en dirección radial
T=np.zeros(len(r))                  #Temperatura [K]
for i in np.arange(len(r)):
  if r[i]<Rn:
    T[i]=Ts+(1./(3.*km))*(am/2.*(R**2-Rn**2)
    +(an-am)*Rn**3*(1./Rn-1./R))+(an/(6.*kn))*(Rn**2-r[i]**2)
    plt.plot(r[i]*1e-3,T[i],'k.')
  else:
    T[i]=Ts+1./(3.*km)*(am/2*(R**2-r[i]**2)+(an-am)*Rn**3*(1./r[i]-1./R))
    plt.plot(r[i]*1e-3,T[i],'r.')

plt.xlabel('Radio terrestre [km]')
plt.ylabel('T(r) [K]')

plt.show()

#----------------------------------------------------------------------
#Caso para un planeta homogeneo y sin estratificación
#----------------------------------------------------------------------

#Número de partículas - Energía de cada isotopo y total - Fuente de calor

N=np.zeros(len(masa_n))    #N° de particulas de la tierra actualmente
No=np.zeros(len(masa_n))   #N° de particulas de la tierra hace 4x10^9 años
e=np.zeros(len(masa_n))    #Energia de cada isotopo de la tierra
e_sum=np.zeros(1)          #Energia total de la tierra

for i in np.arange(len(masa_n)):
    N[i]=(masa[i]/As[i])*Na
    No[i]=N[i]*np.exp(lamb[i]*t)
    e[i]=E[i]*lamb[i]*No[i]
    e_sum=e_sum+e[i]

a=e_sum/Vt                      #Fuente de calor de la tierra [J/m^3s]

#Perfil de calor y de temperatura

#Temperatura superficial

qss  = a/3.*R
Tss=(qss/(eps*sigma)+To**4)**(1./4.)

n=100
r=np.arange(0,R,R/n)                #[m] Posiciones en dirección radial
r[0]=1e-10                          #porque en algun momento divido por r

q=np.zeros(len(r))                  #[J/m²s] Calor por unidad de area y tiempo
TT=np.zeros(len(r))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for j in np.arange(len(r)):
    
    if r[j]<=Rn:
        q[j]=a*r[j]/3.
        TT[j]=Tss+a/(6*k)*(R**2-r[j]**2)
        
        ax1.plot(r[j]*1e-3,q[j],'k.')
        ax1.set_xlabel("Radio terrestre [km]")
        ax1.set_ylabel("q(r) [J/m^2s]")
        ax1.ticklabel_format(useOffset=False)
        ax1.set_xlim(0,R*1e-3)
        
        ax2.plot(r[j]*1e-3,TT[j],'k.')
        ax2.set_xlabel("Radio terrestre [km]")
        ax2.set_ylabel("T(r) [K]")
        ax2.ticklabel_format(useOffset=False)
        ax2.set_xlim(0,R*1e-3)
    
    else:
        q[j]=a*r[j]/3.
        TT[j]=Tss+a/(6*k)*(R**2-r[j]**2)
        
        ax1.plot(r[j]*1e-3,q[j],'r.')
        ax2.plot(r[j]*1e-3,TT[j],'r.')

plt.gcf().axes[0].yaxis.get_major_formatter().set_scientific(False)
plt.show()