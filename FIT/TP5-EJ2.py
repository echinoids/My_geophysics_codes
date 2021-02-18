# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 23:14:45 2021

@author: steph
"""
'''
A partir de los datos que se encuentran en el archivo, estimar
el espesor elastico de la litosfera en esta region y el momento
aplicado sobre la placa litosferica que subducta. 
Utilizar los siguientes valores para los calculos E = 50 GPa,
v= 0;25, rho_m = 3300 kg/m^3, rho_w = 1000 kg/m^3 y g = 10 m/s.

'''

#Importo librerias

import numpy as np
import matplotlib.pyplot as plt

dat=np.loadtxt('fosa.txt')

#Constantes
E=50*1e9    #Modulo de young [Pa]
v=0.25      #Modulo de poisson
rho_m=3300  #Densidad del manto [kg/m^3]
rho_w=1000  #Densidad del agua [kg/m^3]
g=10        #Gravedad [m/s^2]

x=dat[:,0]*1e3  #Extensión horizontal [m]
w=dat[:,1]*1e3  #Deflexión [m]

#Generación de la malla y extensiones

w_max=np.max(w)
idx=np.where(w==w_max)
xb=x[idx]

alpha_ñuflo=(4/np.pi)*xb
beta_ñuflo=w_max*(np.exp(np.pi/4)/np.sin(np.pi/4))

N=50
M=50

alpha_ext=np.linspace(0.5*alpha_ñuflo,1.5*alpha_ñuflo,N)
beta_ext=np.linspace(0.5*beta_ñuflo,1.5*beta_ñuflo,M)

#Funcion de costo e iteracion para hallar el aplha y beta en cada caso

phi=np.zeros((N,M))
X=len(x)
w_new=np.zeros(X)

for i in np.arange(N):
    for j in np.arange(M):
        for k in np.arange(X):
            w_new[k]=beta_ext[j]*np.exp(-x[k]/alpha_ext[i])*np.sin(x[k]/alpha_ext[i])
        #end
        phi[i,j]=np.sum((w-w_new)**2)
    #end
#end

ia,ib=np.where(phi==np.min(phi))

print('Los parámetros óptimos son: alpha={} [km], beta={} [km]'
      .format(alpha_ext[ia]*1e-3,beta_ext[ib]*1e-3))

#Valor del espesor elastico y el momento aplicado

D=(alpha_ext[ia]**4)*(rho_m-rho_w)*(g/4)
h=((D/E)*(12*(1-v**2)))**(1/3)
M=-(2*beta_ext[ib]*D)/alpha_ext[ia]**2

print(' ')
print('El espesor elastico es h={} [m] y el momento es Mo={} [Pam^3]'.format(h,M))

#Graficos

#Datos y datos ajustados

fig=plt.figure()
ax=fig.add_subplot(111)    # The big subplot
ax1= fig.add_subplot(211)
ax2=fig.add_subplot(212)

# Remover los ticks y ejes del subplot grande 
ax.spines['top'].set_color('none'); ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none'); ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

w_optimo=np.zeros(len(x))
for i in range(len(x)):
    w_optimo[i]=beta_ext[ib]*np.exp(-x[i]/alpha_ext[ia])*np.sin(x[i]/alpha_ext[ia])
    
ax1.plot(x*1e-3,w*1e-3,'k*'),ax1.plot(xb*1e-3,w_max*1e-3,'r+',label='$w_{max}(x_b)$')
ax1.legend()
ax2.plot(x*1e-3,w*1e-3,'k*');ax2.plot(x*1e-3,w_optimo*1e-3,'-r',label='Ajuste optimo')
ax2.legend()

ax.set_title('Perfil de deflexion y su ajuste')
ax.set_xlabel('Extensión horizontal [km]')
ax.set_ylabel('Deflexion [km]')

plt.show()

#Malla

xx,yy=np.meshgrid(beta_ext,alpha_ext)
z=np.ones_like(xx)*100

plt.figure(figsize=(10,7),facecolor='w')
plt.pcolormesh(xx*1e-3,yy*1e-3,z,edgecolors='w',cmap="seismic",
               linewidth=1,shading='auto')
plt.title('Mallado')
plt.xlabel(r'$\beta$ [km]')
plt.ylabel(r'$\alpha$ [km]')
plt.show()

#Función de costo y parámetros alfa y beta

plt.figure(figsize=(10,7),facecolor='w')
plt.plot(beta_ext[ia]*1e-3,alpha_ext[ib]*1e-3,'k*')
cb=plt.pcolor(xx*1e-3,yy*1e-3,phi*1e-3/1e6,shading='auto',cmap='RdGy')
clb=plt.colorbar(cb)
clb.ax.text(2.5, 0.5, r'$\phi$($\beta$,$\alpha$)', fontsize=10, rotation=-90, va='center') 
clb.set_label('$x10^6$', labelpad=-30, y=1.06, rotation=0)
plt.title('Matriz de función objetivo')
plt.xlabel(r'$\beta$ [km]')
plt.ylabel(r'$\alpha$ [km]')
plt.show()