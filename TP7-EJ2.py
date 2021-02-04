# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:38:39 2021

@author: steph
"""
import numpy as np
import matplotlib.pyplot as plt

#Constantes
E=21.6*1e9  #Modulo de young [Pa]
E1=E
E2=29.4*1e9
ETA=0.14*1e9 #Viscosidad [Pa.s]

#Maxwell - kelvin - zener

tau=ETA/E
tau_e=ETA/E2
tau_sigma=ETA/(E1+E2)
Er=(E1*E2)/(E1+E2)

w=np.arange(1,1e4,100)

E_real_m=np.zeros(len(w))
E_im_m=np.zeros(len(w))

E_real_k=np.zeros(len(w))
E_im_k=np.zeros(len(w))

E_real_z=np.zeros(len(w))
E_im_z=np.zeros(len(w))

#Factores de disipación Q^-1
Q_m=np.zeros(len(w))
Q_k=np.zeros(len(w))
Q_z=np.zeros(len(w))

for i in range(len(w)):
    
    E_real_m[i]=(ETA*tau*w[i]**2)/(1+tau**2*w[i]**2)
    E_im_m[i]=(ETA*w[i])/(1+tau**2*w[i]**2)
    Q_m[i]=E_im_m[i]/E_real_m[i]

    E_real_k[i]=E
    E_im_k[i]=ETA*w[i]
    Q_k[i]=E_im_k[i]/E_real_k[i]

    E_real_z[i]=(Er+Er*tau_e*tau_sigma*w[i]**2)/(1+tau_sigma**2*w[i]**2)
    E_im_z[i]=((tau_e-tau_sigma)*Er*w[i])/(1+tau_sigma**2*w[i]**2)
    Q_z[i]=E_im_z[i]/E_real_z[i]
        
real=[E_real_m,E_real_k,E_real_z]
im=[E_im_m,E_im_k,E_im_z]
name=['Ensayo din. (Modelo Maxwell)','Ensayo din. (Modelo Kelvin)',
      'Ensayo din. (Modelo Zener)']

fig, axes = plt.subplots(1, len(real),figsize=(10,5))
fig.subplots_adjust(hspace=1, wspace=0.4)

for ax,real,im,name in zip(axes.flat,real,im,name):
    ax.plot(w,real,'k')
    ax.plot(w,im,'r')
    ax.set_title(name)
    ax.set_xlabel("$\omega$")
    ax.set_ylabel('Amp []')
    ax.ticklabel_format(useOffset=False)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(1,1e4)
    
plt.show()

Q=[Q_m,Q_k,Q_z]
name=['Q^-1 (Modelo Maxwell)','Q^-1 (Modelo Kelvin)',
      'Q^-1 (Modelo Zener)']

fig, axes = plt.subplots(1, len(Q),figsize=(10,5))
fig.subplots_adjust(hspace=1, wspace=0.4)

for ax,Q,name in zip(axes.flat,Q,name):
    ax.plot(w,Q,'k')
    ax.set_title(name)
    ax.set_xlabel("$\omega$")
    ax.set_ylabel('Amp []')
    ax.ticklabel_format(useOffset=False)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(1,1e4)
    
plt.show()