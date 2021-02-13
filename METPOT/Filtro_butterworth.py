#Importo módulos
import numpy as np
import matplotlib.pyplot as plt
#
# #Defino el número de onda en cada dirección y el número de onda radial
# def kr (nx,ny,delx,dely):
#     #
#     kx = np.zeros(nx, dtype=float)
#     #
#     #Evaluando cada k y preguntando si k >= Nyquist == pi
#     #
#     dk = (2 * np.pi / nx)
#     #
#     for j in np.arange(nx):
#         kx[j] = dk * j
#         if kx[j] >= np.pi:
#             kx[j] = kx[j] - 2.0 * np.pi
#     #
#     # Para mayor facilidad, para calcular ky, definiremos:
#     #
#     ky,ny = kx,nx
#     #
#     #Calculo el numero de onda radial
#     kr = np.zeros(nx*ny).reshape(nx,ny)
#     #
#     for i in np.arange(nx):
#         kr[i, :] = np.sqrt((kx[i]/delx)** 2 + (ky[:]/dely)** 2)
#     #
#     return kr
#
#Defino una subrutina para el filtro de Butterworth en el dominio de las frecuencias
#
def butter (nx,ny,n,kr,kc):
    BW=np.zeros(nx*nx).reshape(nx,ny)
    for j in range (nx):
        for i in range (ny):
            BW[i,j] = 1/np.sqrt(((1+(kr[i,j]/kc)**n)))
    return BW
#Introducimos los datos sintéticos
#
T = np.loadtxt("http://carina.fcaglp.unlp.edu.ar/~jgomez/academic/mpp/notebooks/sintetico_pole_reduction.dat").T # x [m], y [m], T [nT]
#
#Defino parámetros de longiud y espaciamiento en x e y, a fin de que sen matrices 100x100
L      = 7000         #m
nx, ny = (100,100)    #m
#
#Busco un espaciamiento que me permita ordenar los datos en dicha matriz
dx     = 2*L/nx
dy     = 2*L/ny
#
#Defino vector de muestras equiespaciadas desde -L a L, con espaciamiento nx y ny
xx    = np.linspace(-L, L, nx)  #m
yy    = np.linspace(-L, L, ny)  #m
#
#Defino la grilla
x, y  = np.meshgrid(xx, yy)
#
dato  = np.reshape(T,(ny,nx))
#
#Defino el ruido aleatorio con distribución Normal estándar de media 0 y deviación 1 para sumarlo a la señal
#multiplico por un factor de amplitud que será un porcentaje de la amplitud maxima de la anomalía
#
#Porcentaje de ruido
noise_p= 2
#
#Ruido
noise = (np.max(dato))*(noise_p/100)*np.random.normal(0,1,len(T))
Noise = np.reshape(noise,(nx,ny))
#
#Añado ruido
noisydat = dato + Noise
#
#Voy al dominio transformado
DATO = np.fft.fft2(noisydat)
#
#Defino deltas de longitud en cada dirección delx y dely
delx  =xx[2]-xx[1]
dely  =yy[2]-yy[1]
#
#Defino el número de onda de otra manera mas cómoda
#Puede ser por la subrutina
#kr    = kr(nx,ny,delx,dely)
#
#O por esta manera un poco mas rápida
kx = (2*np.pi)*np.fft.fftfreq(nx)#/delx #Le doy dimensiones al número de onda
ky = (2*np.pi)*np.fft.fftfreq(ny)#/dely

kr = np.zeros(nx*ny).reshape(nx,ny)
for i in np.arange(nx):
    kr[i,:] = np.sqrt(kx[i]**2+ky[:]**2)
#
#Defino los parámetros del filtro n "orden" y Kc "número de onda de corte"
n=8
kc=np.pi/400
filtro = butter(nx,ny,n,kr,kc)
#
#Aplicación del filtro de Butterworth
DATO_fil = DATO*filtro
#
#Vuelvo al dominio espacial
dato_fil = np.real(np.fft.ifft2(DATO_fil))
#
#Gráfico de los datos anteriores
#------------------------------------------------Figura 1---------------------------------------------
# plt.figure(1,figsize=(5,5))
# plt.title(r'$\Delta T$ del dipolo')
# plt.xlabel("$m$")
# plt.ylabel("$m$")
# #
# plt.imshow(dato,origin="lower",cmap="RdYlBu_r",extent=[-L,L,-L,L],interpolation='nearest')
# plt.colorbar(orientation='horizontal',pad=0.2).set_label('[nT]', labelpad=5, y=0.5, rotation=0)
# plt.contour(xx,yy,dato,15,linewidths=0.5,colors='k')
# #
# plt.show()
#--------------------------------------------------Figura 2------------------------------------------
# plt.figure(2,figsize=(5,5))
# plt.title(r'$\Delta T$ del dipolo')
# plt.xlabel("$m$")
# plt.ylabel("$m$")
# #
# plt.imshow(noisydat,origin="lower",cmap="RdYlBu_r",extent=[-L,L,-L,L],interpolation='nearest')
# plt.colorbar(orientation='horizontal',pad=0.2).set_label('[nT]', labelpad=5, y=0.5, rotation=0)
# plt.contour(xx,yy,noisydat,10,linewidths=0.5,colors='k')
# #
# plt.show()
#------------------------------------------------Figura 1---------------------------------------------
plt.figure(3,figsize=(5,5))
plt.title("$|Filtro\,\,Butterworth|$")
plt.xlabel("$k_x$")
plt.ylabel("$k_y$")
#
plt.xticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
plt.yticks([0,np.pi,2*np.pi],[r'$0$',r'$\pi$',r'$2\pi$'])
extent=[0,2*np.pi,0,2*np.pi]
#
plt.imshow(np.abs(np.fft.fftshift(filtro)),origin="lower",extent= extent,cmap="gray",interpolation='none')
plt.colorbar(orientation='horizontal',pad=0.2)
plt.show()
#------------------------------------------------Figura 3---------------------------------------------
# plt.figure(4,figsize=(5,5))
# plt.title(r'$\Delta T$ del dipolo')
# plt.xlabel("$m$")
# plt.ylabel("$m$")
# #
# plt.imshow(dato_fil,origin="lower",cmap="RdYlBu_r",extent=[-L,L,-L,L],interpolation='nearest')
# plt.colorbar(orientation='horizontal',pad=0.2).set_label('[nT]', labelpad=5, y=0.5, rotation=0)
# plt.contour(xx,yy,dato_fil,10,linewidths=0.5,colors='k')
# #
# plt.show()
# #------------------------------------------------Figura 3---------------------------------------------
# plt.figure(5,figsize=(5,5))
# plt.title(r'$\Delta T$ del dipolo')
# plt.xlabel("$m$")
# plt.ylabel("$m$")
# #
# plt.imshow(dato-dato_fil,origin="lower",cmap="RdYlBu_r",interpolation='nearest',extent=[-L,L,-L,L],vmin=np.amin(dato), vmax=np.amax(dato))
# plt.colorbar(orientation='horizontal',pad=0.2).set_label('[nT]', labelpad=5, y=0.5, rotation=0)
# plt.contour(xx,yy,dato-dato_fil,10,linewidths=0.5,colors='gray')
# #
# plt.show()

dat = [dato,dato_fil,noisydat,dato-dato_fil]
names   = [r'$\Delta T$ del dipolo [nT]',r'$\Delta T$ filtrado [nT]',r'$\Delta T$ ruidoso [nT]',r'$\Delta T$ residual [nT]']
#
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
fig.subplots_adjust(hspace=0.1, wspace=0)
#
for ax, dat,names in zip(axes.flat, dat,names):
    im=ax.imshow(dat, origin="lower",interpolation='nearest',extent=[-L,L,-L,L],vmin=np.amin(dato), vmax=np.amax(dato), cmap='RdYlBu_r')
    ax.set_title(names)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.contour(xx, yy,dat, 10, linewidths=0.5, colors='gray')
    fig.colorbar(im,ax=ax, orientation='horizontal',pad=0.2,shrink=0.335,spacing='uniform')
#
plt.show()