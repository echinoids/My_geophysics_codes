# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:29:55 2021

@author: steph
"""
import numpy as np
import matplotlib.pyplot as plt

Em=2860e3                        #Espesor del manto [m]
Ec=30e3                          #Espesor de la corteza [m]
Rn=3480e3                        #Radio del núcleo [m]
Rm=(Rn+Em)                       #Radio del manto [m]
Rt=(Rm+Ec)                       #Radio de la Tierra [m]
Na=6.022e23                      #Número de avogadro [1/mol]
Vn= (4./3.)*np.pi*Rn**3          #Volumen del núcleo [m^3]
Vm= (4./3.)*np.pi*(Rm**3-Rn**3)  #Volumen del manto [m^3]
Vc= (4./3.)*np.pi*(Rt**3-Rm**3)  #Volumen de la corteza [m^3]
E=76.245e-13                     #Energía [Joules]
lamb=9.85e-10                    #Constante de decaimiento [1/año]
masa_tierra=5.97e24              #Masa total de la Tierra [kg]
As=0.230                         #Masa atómica del elemento [kg]
yrtosec=365*24*60*60             #Ctte para llevar de año a segundos
lamb=lamb/yrtosec                #Constante de decaimiento [1/s]

#----------------------------------------------------------------------
#Defino funcion de calor especificamente para la superficie terrestre
#----------------------------------------------------------------------

def q_c(an,am,ac,Rn,Rm,Rt):
    
    '''
    Del inciso 1.a tomamos la expresión para la función
    que describe el flujo de calor en el intervalo Rm < r < Rt
    ya que nos indican en valor que debe tener para r = Rt
    
    '''
    
    qc= (1/3)*(ac*Rt+(1/Rt)**2*((Rm**3)*(am-ac)+(Rn**3)*(an-am)))
    return qc

#----------------------------------------------------------------------
#Casos para un planeta con estratificación (nucleo-manto-corteza)
#----------------------------------------------------------------------

casos=3
prop=[1.16e-5,1.16e-8,1.16e-11]  #Proporcion de la masa total de la Tierra
Pn=[30,0,25]                     #Porcentajes de la masa total que ocupa
Pm=[70,70,50]                    #el elemento entre el nucleo (Pn), manto (Pm)
Pc=[0,30,25]                     #y corteza (Pc)
an=np.zeros(casos)
am=np.zeros(casos)               #Vectores de fuentes de calor para cada zona
ac=np.zeros(casos)

for i in range(casos):

    masa_n=prop[i]*masa_tierra*(Pn[i]/100)     #Masa del nucleo [kg]
    masa_m=prop[i]*masa_tierra*(Pm[i]/100)     #Masa del manto [kg]
    masa_c=prop[i]*masa_tierra*(Pc[i]/100)     #Masa de la corteza [kg]

    Nn=(masa_n/As)*Na           #Numero de particulas del nucleo
    Nm=(masa_m/As)*Na           #Numero de particulas del manto
    Nc=(masa_c/As)*Na           #Numero de particulas de la corteza

    en=E*lamb*Nn                #Energia del nucleo
    em=E*lamb*Nm                #Energia del manto
    ec=E*lamb*Nc                #Energia de la corteza

    an[i]=en/Vn                 #Fuente de calor del nucleo[J/m^3s]
    am[i]=em/Vm                 #Fuente de calor del manto [J/m^3s]
    ac[i]=ec/Vc                 #Fuente de calor de la corteza [J/m^3s]

    #Calculo el calor en la superficie de la Tierra

    qc=q_c(an[i], am[i], ac[i], Rn, Rm, Rt)
    
    print('')
    print('Caso N° {}'.format(i+1))
    print('')
    print('El valor del flujo de calor en la superficie es q_c(Rt)={:e} [J/m^2s]'.format(qc))
    print('')

#----------------------------------------------------------------------
#Visualización del flujo en la corteza  Rm < r  < Rt
#----------------------------------------------------------------------

#Tomando los datos del caso correcto (Caso 2)

a_n=an[1];a_m=am[1];a_c=ac[1]    #Fuentes de calor
r=np.arange(Rm,Rt+100,100)       #[m] Posiciones en dirección radial
q=np.zeros(len(r))               #[J/m²s] Calor por unidad de area y tiempo

fig1, ax1 = plt.subplots()

for i in np.arange(len(r)):
    
    q[i]=(1/3)*(a_c*r[i]+(1/r[i])**2*((Rm**3)*(a_m-a_c)+(Rn**3)*(a_n-a_m)))
    
    ax1.plot(r[i]*1e-3,q[i],'k.')
    ax1.set_xlabel("Radio terrestre [km]")
    ax1.set_ylabel("$q_c(r)$ [J/m^2s]")
    ax1.set_xlim(Rm*1e-3,Rt*1e-3)

plt.show()