import numpy as np
import matplotlib.pyplot as plt

def convolve (x,y):
    '''
    Genero convolución inmediata
    '''
    
    xconv=np.convolve(x,y,mode='same')
    
    return xconv


def ricker(cfreq,phase,dt,wvlt_length):
    '''
    Calculate a Ricker wavelet
    
    Usage:
    --------
    t, wvlt = wvlt_ricker(cfreq, phase, dt, wvlt_length)
    
    cfreq:central frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    '''

    t_max = wvlt_length*0.5
    t_min = -t_max
    t = np.arange(t_min,t_max,dt)
    t = np.linspace(-wvlt_length/2.0, (wvlt_length-dt)/2.0, np.int(wvlt_length/dt))
    wvlt = (1.0 - 2.0*(np.pi**2.0)*(cfreq**2.0)*(t**2.0))*np.exp(-(np.pi**2.0)*(cfreq**2.0)*(t**2.0))
        
    return t, wvlt

#---------------------------------------------------------------------#

dat=np.loadtxt('datos.txt')
dat=np.array(dat)

TWT=dat[:,0]
Z=dat[:,1]
CR=dat[:,2]
RHO=dat[:,3]


wvlt_length=0.5      #Wavelet length in seconds
wvlt_phase=0.0       #Wavelet phase in degrees
wvlt_cfreq=30.0      #Ricker wavelet central frequency
dt=0.004             #Intervalo de muestreo

wvlt_t, wvlt_amp = ricker(wvlt_cfreq,wvlt_phase,dt,wvlt_length)

synth=convolve(CR,wvlt_amp)

# =================================================================
# Gráficos
# =================================================================


fig = plt.figure(figsize=(20, 10))

ax = fig.add_subplot(2, 1, 2)
ax.plot(TWT,Z)
ax.set_xlim(np.min(TWT), np.max(TWT))
ax.set_ylim(np.min(Z),np.max(Z))
ax.set_title('Impedancias vs TWT')
ax.set_ylabel("Amplitud")
ax.set_xlabel("TWT")

ax = fig.add_subplot(2, 2, 1)
ax.plot(TWT,CR)
ax.set_xlim(np.min(TWT), np.max(TWT))
ax.set_ylim(np.min(CR),np.max(CR))
ax.set_title('Coeficientes de reflexión vs TWT')
ax.set_ylabel("Amplitud")
ax.set_xlabel("TWT")

ax = fig.add_subplot(2, 2, 2)
ax.plot(TWT,RHO)
ax.set_xlim(np.min(TWT), np.max(TWT))
ax.set_ylim(np.min(RHO),np.max(RHO))
ax.set_title('Densidad vs TWT')
ax.set_ylabel("Amplitud")
ax.set_xlabel("TWT")

plt.plot()

fig = plt.figure(figsize=(20, 10))

ax = fig.add_subplot(1, 2, 1)
ax.plot(synth,TWT)
ax.xaxis.tick_top()
ax.set_ylim(np.min(TWT), np.max(TWT))
plt.gca().invert_yaxis()
ax.set_xlim(np.min(synth)-0.1,np.max(synth)+0.1)
ax.set_title('Sismograma sintético')
ax.set_ylabel("TWT")
ax.set_xlabel("Amp")

ax = fig.add_subplot(1, 2, 2)
ax.plot(wvlt_amp,wvlt_t)
ax.xaxis.tick_top()
ax.set_ylim(np.min(wvlt_t), np.max(wvlt_t))
ax.set_xlim(np.min(wvlt_amp),np.max(wvlt_amp))
ax.set_title('Ondicula de Ricker')
ax.set_ylabel("t")
ax.set_xlabel("Amp")

plt.plot()