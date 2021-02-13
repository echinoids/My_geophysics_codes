#Importo módulos
import matplotlib.pyplot as plt
import numpy as np
import numpy
import sys
from scipy.interpolate import griddata
numpy.set_printoptions(threshold=sys.maxsize)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#
#Rutina de anomalía de esfera
#
def g_vert(x,y,z,X0,R,drho):
  #Constante de gravitación universal
  G=6.67408e-11  #m^3/kg.s^2
  #
  #Constante de conversión de m/s^2 a mGal
  conv=1e5
  #
  #Amplitud de la contribución
  ctte=((4/3)*np.pi*(R**3)*drho*G)*conv
  #
  #Contribución
  x0,y0,z0=X0
  g_z=ctte*((z-z0)/(np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3))
#
  return g_z
#
def kr(nx):
    #
    kx = np.zeros(nx, dtype=float)
    #
    #Evaluando cada k y preguntando si k >= Nyquist == pi
    #
    dk = (2 * np.pi / nx)
    #
    for j in np.arange(nx):
        kx[j] = dk * j
        if kx[j] >= np.pi:
            kx[j] = kx[j] - 2.0 * np.pi
    #
    # Para mayor facilidad, para calcular ky, definiremos:
    #
    ky, ny = kx, nx
    #
    #Calculo el numero de onda radial
    kr = np.zeros(nx*ny).reshape(nx,ny)
    #
    for i in np.arange(nx):
        kr[i, :] = np.sqrt(kx[i] ** 2 + ky[:] ** 2)
    #
    return kr
    #
#Defino función de continuación ascentente o descendente (Continuación analitica)
# Z > 0 para continuación ascendente
# Z < 0 para continuación descendente
#
def cont_ac(kr,z):
  #
  H = np.exp(-kr*z)
  #
  return H
#
#Defino parámetros de longiud y espaciamiento en x e y:
L      = 1100    #m
nx, ny = (81,81) #m
#
#Defino vector de muestras equiespaciadas desde -L a L, con
#espaciamiento nx y ny:
xx     = np.linspace(-L, L, nx)  #m
yy     = np.linspace(-L, L, ny)  #m
#
#Defino la grilla:
x, y   = np.meshgrid(xx, yy)
#
#Defino arreglo de ceros con la misma longitud que x:
z      = np.zeros_like(x)
#
#Parámetros
#
#Contraste de densidad
drho  =500              #kg/m^3
#Radio de la esfera
R     =200              #m
#Ubicación del/los cuerpos (x0,y0,z0)
X01   =(240,0,-400)     #m
X02   =(-240,0,-400)    #m
#
#Contribuciones
g_z1  =g_vert(x,y,z,X01,R,drho)  #mGal
g_z2  =g_vert(x,y,z,X02,R,drho)  #mGal
g_zt  =g_z1+g_z2                 #mGal
#
#Aplico la continuación analítica a los datos
#
#Defino el filtro en frecuencia
z  = -100
kr = kr(nx)
filtro = np.fft.fftshift(cont_ac(kr,z))
#
#Transformo mis datos del dominio espacial al dominio de los números de onda
G_ZT      = np.fft.fft2(g_zt)
#
#Centralizo
G_ZT2 = np.fft.fftshift(G_ZT)
#
#Aplico el filtro en frecuencia
G_ZT_freq = filtro*G_ZT2
#
#Descentralizo
G_ZT_freq2 = np.fft.ifftshift(G_ZT2)
#
#Antitransformo para ir al dominio espacial
g_zt_fil  = np.real(np.fft.ifft2(G_ZT_freq2))

#Grafico los datos anteriores
#------------------------------------------------Figura 1---------------------------------------------
plt.figure(1,figsize=(5,5))
plt.title(r'$\Delta g_z$ de ambas esferas')
plt.xlabel("$m$")
plt.ylabel("$m$")
plt.imshow(g_zt,origin="lower",cmap="coolwarm",interpolation='nearest')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.contour(g_zt_fil,10,linewidths=0.5,colors='k')
plt.show()
#--------------------------------------------------Figura 2------------------------------------------
plt.figure(2,figsize=(5,5))
plt.title("$|k_r|$")
plt.xlabel("$k_x$")
plt.ylabel("$k_y$")
#
plt.xticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
plt.yticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
extent=[0,2*np.pi,0,2*np.pi]
#
plt.imshow(kr,origin="lower",extent= extent,cmap="gray",interpolation='none')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.show()
#-----------------------------------------------Figura 3---------------------------------------------------
plt.figure(3,figsize=(5,5))
plt.title("$Filtro$")
plt.xlabel("$k_x$")
plt.ylabel("$k_y$")
#
plt.xticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
plt.yticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
extent=[0,2*np.pi,0,2*np.pi]
plt.imshow(np.abs(filtro),origin="lower",extent=extent,cmap="gray",interpolation='none')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.show()
#-----------------------------------------Figura 4---------------------------------------------------
plt.figure(4,figsize=(5,5))
plt.title("$G_{ZT}$")
plt.xlabel("$k_x$")
plt.ylabel("$k_y$")
#
plt.xticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
plt.yticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
extent=[0,2*np.pi,0,2*np.pi]
#
plt.imshow(np.abs(G_ZT),origin="lower",extent=extent,cmap="gray",interpolation='none')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.show()
#------------------------------------Figura 5-------------------------------------------------
plt.figure(5,figsize=(5,5))
plt.title("$G_{ZT}$ filtrada")
plt.xlabel("$k_x$")
plt.ylabel("$k_y$")
#
plt.xticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
plt.yticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
extent=[0,2*np.pi,0,2*np.pi]
#
plt.imshow(np.abs(G_ZT_freq),origin="lower",extent=extent,cmap="gray",interpolation='none')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.show()
#-------------------------------------------Figura 6-------------------------------------------------
plt.figure(6,figsize=(5,5))
plt.title("$g_{zt}$ filtrada")
plt.xlabel("$m$")
plt.ylabel("$m$")
#
plt.imshow(np.abs(g_zt_fil),origin="lower",cmap="coolwarm",interpolation='nearest')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.contour(np.abs(g_zt_fil),10,linewidths=0.5,colors='k')

plt.show()
