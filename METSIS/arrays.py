import numpy as np
import matplotlib.pyplot as plt
import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)
np.seterr(divide='ignore', invalid='ignore')


nr = [12,24,12,6,6,4,12]  # nro de receptores
dx = [4.0,2.0,2.0,2.0,4.0,12.5,4.17]  # distancia entre receptores
nk = 1500  # puntos del grafico

# frecuencia de Nyquist espacial
fnyq=np.zeros(len(dx))
for i in range (len(dx)):
    fnyq[i]=1.0/(2.0*dx[i])

#mutliplicar por un flotante genera errores, por ello se usan este tipo de
#escrituras en lista para efectuar mutliplicaciones de constantes reales por arrays

# nro de onda maximo
kmax = [i * (2.0) for i in fnyq]

def k (kmax,dx,nk):
    k= dx* np.linspace(0, kmax, nk)  # rango de valores de k
    return k

K=np.zeros(shape=(nk,len(dx)))
for i in range (len(dx)):
    kw= k(kmax[i],dx[i],nk)
    K[:,i]=kw

R=np.zeros(shape=(nk,len(dx)))
for i in range (len(nr)):
    R[:,i] = np.sin(nr[i] * np.pi * K[:,i]) / (nr[i] * np.sin(np.pi * K[:,i]))  # respuesta
    R[np.isnan(R)] = 1.0
#
fig = plt.figure()
fig.subplots_adjust(hspace = 1.5, wspace=0.5)
#
for i in range (len(dx)):
    ax = fig.add_subplot(3, 3, i+1)
    ax.set_xlim([0, dx[i] * kmax[i]])
    ax.set_ylim([0, 1.0])
    ax.set_xlabel(r'dx.k=dx/$\lambda$')
    ax.set_ylabel(r'R')
    ax.set_title('Arreglo {}'.format(i+1))
    ax.plot(K[:,i], np.abs(R[:,i]))
#
plt.show()
