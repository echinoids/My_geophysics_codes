# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 12:25:02 2020

@author: RY11675
"""

def ConveccionTierraCapaLimite(R,Rc,T0,Rac,C_c,C_m,n0,fm,fc,t0,tf,dt,Tc0,Tm0):
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
    import numpy as np
    # Ctes
    c2 = 2**(4/3)/Rac**(1/3)
    
    # datos Nucleo
    ka_c = C_c[0]
    k_c = C_c[1]
    alfa_c = C_c[2]
    nu_c = C_c[3]
    g_c = C_c[4]
    Vc = 4*np.pi*Rc**3/3 
    Ac = c2*k_c*(alfa_c*g_c/ka_c/nu_c)**(1/3)
    
    # datos Manto
    ka_m = C_m[0]
    k_m = C_m[1]
    alfa_m = C_m[2]
    nu_m = C_m[3]
    g_m = C_m[4]
    Vm = 4*np.pi*(R**3-Rc**3)/3
    Am = c2*k_m*(alfa_m*g_m/ka_m/nu_m)**(1/3)
    
    # Defino ctes para los isotopos
    eV = 1.6e-19
    E = np.array([47.7, 43.9, 40.5, 0.71])*eV*1.e6 # Joules
    y = 365*24*3600
    l = np.array([1.55e-10, 9.85e-10, 4.94e-11, 5.54e-10])/y # 1/seg
    
    # Calculo concentración de isotopos al tiempo 0 en el manto y el nucleo
    n0_m = fm*n0
    n0_c = fc*n0
    
    # defino el vector de tiempo
    t = np.arange(t0,tf+dt,dt)
    # inicializo arreglos
    Tc = np.zeros(len(t))
    Tm = np.zeros(len(t))
    
    # Condicione iniciales
    Tc[0]=Tc0
    Tm[0]=Tm0
    
    for k,tk in enumerate(t[:-1]):
        # calculo la fuente de calor radiogenico
        n_m = n0_m*np.exp(-l*tk*y)
        am = np.sum((n_m*l)*E)/Vm 
        n_c = n0_c*np.exp(-l*tk*y)
        ac = np.sum((n_c*l)*E)/Vc
        # Calculo temperatura en la interface Nucleo-Manto
        Tmc = (Ac**(3/4)*Tc[k]+Am**(3/4)*Tm[k])/(Ac**(3/4)+Am**(3/4))
        # Calculo las temperatura del manto y nucleo usando Euler
        AAc = 3*ka_c*Ac/k_c/Rc
        dTc = ka_c*ac/k_c-AAc*(Tc[k]-Tmc)**(4/3) # K/seg
        Tc[k+1] = Tc[k] + dt * dTc * y # K
        AAm = Am*3*ka_m/k_m/(R**3-Rc**3)
        dTm = ka_m*am/k_m + Rc**2*AAm*(Tmc-Tm[k])**(4/3)-R**2*AAm*(Tm[k]-T0)**(4/3)
        Tm[k+1] = Tm[k] + dt * dTm * y
        
    return t,Tc,Tm

