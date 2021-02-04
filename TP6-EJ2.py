# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 20:39:25 2021

@author: steph
"""
##############################################################################
#  PARTE INICIAL
##############################################################################

#Importo modulos

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString

#Constantes para la primera parte
gamma= 1        #Parámetro de Grüneisen
rho=5000        #Densidad [kg/m^3]
kt= 500*1e9     #Incompresibilidad isotérmica [GPa]
Tf_R=1500       #Temperatura para Tf(P=0) [K]
R=6370*1e3      #Radio terrestre [m]
Rc=3480*1e3     #Radio del núcleo [m]
G=6.67430*1e-11 #Constante de gravitación universal [Nm^2/kg^2]

n=100
r=np.arange(0,R,R/n)                 #[m] Posiciones en dirección radial
T0=np.zeros(len(r))                  #Temperatura [K]

fig1, ax = plt.subplots()

for j in np.arange(len(r)):
    
    if r[j]<=Rc:
        T0[j]=Tf_R*np.exp((4./3.)*np.pi*rho**2*G*(gamma-1./3.)/kt*(R**2-r[j]**2))
        
        ax.plot(r[j]*1e-3,T0[j],'k.')
        
    else:
        T0[j]=Tf_R*np.exp((4./3.)*np.pi*rho**2*G*(gamma-1./3.)/kt*(R**2-r[j]**2))
    
        ax.plot(r[j]*1e-3,T0[j],'r.')
        
    ax.set_xlabel("Radio terrestre [km]")
    ax.set_ylabel("Tf(r) [K]")
    ax.ticklabel_format(useOffset=False)
    ax.set_xlim(0,R*1e-3)
    ax.set_title("Temperatura de fusión")
    
plt.show()

##############################################################################
#  PARTE COMPARATIVA
##############################################################################

dat=np.loadtxt('tabla.txt')

#Constantes para la segunda parte

R=6370*1e3                         #Radio terrestre [m]
Na=6.02*1e23                       #Número de avogadro
Vn= (4./3.)*np.pi*Rc**3            #Volumen del núcleo [m^3]
Vm= (4./3.)*np.pi*(R**3-Rc**3)     #Volumen del manto [m^3]
Vt= (4./3.)*np.pi*R**3             #Volumen de la tierra [m^3]
To=270                             #Temperatura inicial [K]
kn=36                              #Constante de resistividad del nucleo [W/mK]
km=6                               #Constante de resistividad del manto [W/mK]
k=8                                #Constante de resistividad de la tierra [W/mK]
sigma=5.67*1e-8                    #Constante de Boltzmann [W/m^2K^4]
eps=1                              #Constante adimensional epsilon
t0=4*1e9                           #Tiempo inicial
t=t0*365*24*3600                   #Tiempo en el pasado [s]

#----------------------------------------------------------------------
#Caso para un planeta con estratificación
#----------------------------------------------------------------------

#Defino vectores y unidades m.k.s

E=dat[:,0]        #Energía [MeV]
lamb=dat[:,1]     #Constante de decaimiento[1/s]
masa=dat[:,2]     #Masa total [kg]
As=dat[:,3]       #Masa atómica de los elementos [kg]

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

#Perfil de temperatura

Ts=(am*R/3.*sigma*eps+((an-am)*Rc**3)/(3.*sigma*eps*R**2)+To**4)**(1/4)

T1=np.zeros(len(r))                  #Temperatura [K]
for i in np.arange(len(r)):
  if r[i]<Rc:
    T1[i]=Ts+(1./(3.*km))*(am/2.*(R**2-Rc**2)
    +(an-am)*Rc**3*(1./Rc-1./R))+(an/(6.*kn))*(Rc**2-r[i]**2)
  else:
    T1[i]=Ts+1./(3.*km)*(am/2*(R**2-r[i]**2)+(an-am)*Rc**3*(1./r[i]-1./R))

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

#Perfil de temperatura

#Temperatura superficial

qss  = a/3.*R
Tss=(qss/(eps*sigma)+To**4)**(1./4.)
TT=np.zeros(len(r))

for j in np.arange(len(r)):
    TT[j]=Tss+a/(6*k)*(R**2-r[j]**2)

#Grafico comparativo

n=100
r=np.arange(0,R,R/n)                #[m] Posiciones en dirección radial
r[0]=1e-10                          #porque en algun momento divido por r

datos=[T1,TT]
name=['Temp. tierra estratificada vs Tf','Temp. tierra no estratificada vs Tf']

fig, axes = plt.subplots(1, len(datos),figsize=(10,5))
fig.subplots_adjust(hspace=1, wspace=0.4)

for ax,datos,name in zip(axes.flat,datos,name):
    ax.plot(r*1e-3,datos,'k')
    ax.plot(r*1e-3,T0,'b')
    ax.set_title(name)
    ax.set_xlabel("Radio [Km]")
    ax.set_ylabel("Temp [K]")
    ax.ticklabel_format(useOffset=False)
    ax.set_ylim(0,np.max(datos))
    ax.set_xlim(0,R*1e-3)

plt.gcf().axes[0].yaxis.get_major_formatter().set_scientific(False)
plt.show()

#Intersección de los graficos anteriores

#Forma 1
first_line = LineString(np.column_stack((r*1e-3, T0)))
second_line = LineString(np.column_stack((r*1e-3, T1)))
intersection = first_line.intersection(second_line)

#Forma 2
index=np.where(T1<=T0)
min_index=np.min(index)
T_min=np.round(T1[min_index],3)
R_min=(R-r[min_index])*1e-3

print(f"La intersección se produce en {R_min} km {T_min} K")