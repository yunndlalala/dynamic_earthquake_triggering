"""
@version:
author:yunnaidan
@time: 2019/04/21
@file: welch_test2.py
@function:
"""
import numpy as np
from scipy import signal
from scipy.integrate import simps
import matplotlib.pyplot as plt
#Generate a test signal, a 2 Vrms sine wave at 1234 Hz, corrupted by 0.001 V**2/Hz of white noise sampled at 10 kHz.

fs = 1e4   #sampling rate, dt = 1/fs
N = 1e5
amp = 2*np.sqrt(2)
freq = 1234.0
noise_power=0.001 * fs / 2
time = np.arange(N) / fs
x = amp*np.sin(2*np.pi*freq*time)
x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)


# Estimate PSD `S_xx_welch` at discrete frequencies `f_welch`
f_welch, S_xx_welch = signal.welch(x, fs=fs,nperseg=1024)
#f_welch_spec, S_xx_welch_spec = signal.welch(x, fs=fs,nperseg=1024,scaling='spectrum')
# Integrate PSD over spectral bandwidth
# to obtain signal power `P_welch`
df_welch = f_welch[1] - f_welch[0]
P_welch = np.sum(S_xx_welch) * df_welch
#P_welch_spec=np.sum(S_xx_welch_spec)* df_welch

###########################
# Compute DFT
Xk = np.fft.fft(x)

# Compute corresponding frequencies
dt = time[1] - time[0]
f_fft = np.fft.fftfreq(len(x), d=dt)

# Estimate PSD `S_xx_fft` at discrete frequencies `f_fft`
T = time[-1] - time[0]
S_xx_fft = ((np.abs(Xk) * dt) ** 2) / T

# Integrate PSD over spectral bandwidth to obtain signal power `P_fft`
df_fft = f_fft[1] - f_fft[0]
P_fft = np.sum(S_xx_fft) * df_fft


#############################
P_exp = (amp / np.sqrt(2)) ** 2
P_exp += noise_power

#############################
E=np.sum(x**2)

print (P_welch)
#print(P_welch_spec)
print (P_fft)
print(P_exp)
print (np.sum(x**2))
print (np.sum(np.abs(Xk) ** 2)/N)
print (np.sum(x**2)*dt/T)
print ('Finish')