# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:30:54 2021

@author: steph
"""
import numpy as np
import matplotlib.pyplot as plt

#Considerar un material con un modulo elastico E = 200 MPa y una viscosidad 
#25 GPahr es comprimido con una tension igual 10 MPa durante 100 horas. 
#Calcular la deformacion luego de 50, 150 y 300 horas utilizando los 
#modelos de Maxwell y de Kelvin-Voigt.

#Constantes
E=200*1e6  #Modulo de young [Pa]
ETA=25*1e9 #Viscosidad [Pa.hora]
sigma0=10*1e6 #Tensión [Pa]
to = 100  #Tiempo para el cual mantengo la tensión constante

#Maxwell y kelvin

t0=int(input('Ingrese el tiempo para calcular la deformación :'))

tau=ETA/E
t=np.arange(0,t0,1)
e_maxwell=np.zeros(len(t))
e_kelvin=np.zeros(len(t))
sigma=np.zeros(len(t))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

for i in range(len(t)):
    
    if t[i]<=to:
        
        sigma[i]=sigma0
        e_maxwell[i]=(sigma0/E)*(1+t[i]/tau)
        e_kelvin[i]=(sigma0/E)*(1-np.exp(-t[i]/tau))
        
        ax1.plot(t[i],e_maxwell[i],'k.')
        ax1.set_title('Deformación para t= {} horas (Modelo de Maxwell)'.format(t0))
        ax1.set_xlabel('Tiempo [horas]')
        ax1.set_ylabel('Deformación []')
        
        ax2.plot(t[i],e_kelvin[i],'k.')
        ax2.set_title('Deformación para t= {} horas (Modelo de Kelvin)'.format(t0))
        ax2.set_xlabel('Tiempo [horas]')
        ax2.set_ylabel('Deformación []')
        
        ax3.plot(t[i],sigma[i],'k.')
        ax3.set_title('Tensión aplicada')
        ax3.set_xlabel('Tiempo [horas]')
        ax3.set_ylabel('Tensión [Pa]')
   
    else:
        
        sigma[i]=0
        e_maxwell[i]=(sigma0*to)/ETA
        e_kelvin[i]=(sigma0/E)*(np.exp((to-t[i])/tau)-np.exp(-t[i]/tau))
        
        ax1.plot(t[i],e_maxwell[i],'g.')
        ax2.plot(t[i],e_kelvin[i],'g.')
        ax3.plot(t[i],sigma[i],'g.')

ax1.plot(t,e_maxwell,'k--',linewidth=0.5)
ax2.plot(t,e_kelvin,'k--',linewidth=0.5)
ax3.plot(t,sigma,'k--',linewidth=0.5)

ax1.set_xlim(0,t0)
ax2.set_xlim(0,t0)
ax3.set_xlim(0,t0)

plt.show()

print(' ')
print('La deformación para t0={} es {} [Pa] con el modelo de Maxwell.'.format(t0,e_maxwell[len(t)-1]))
print(' ')
print('La deformación para t0={} es {} [Pa] con el modelo de Kelvin.'.format(t0,e_kelvin[len(t)-1]))