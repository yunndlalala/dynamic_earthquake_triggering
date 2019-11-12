"""
@version:
author:yunnaidan
@time: 2019/03/27
@file: welch_test.py
@function:
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.integrate import simps

fs = 100
N = 1e4
amp = 2*np.sqrt(2)
freq = 2
noise_power = 0.001 * fs / 2
time = np.arange(N) / fs
x = amp*np.sin(2*np.pi*freq*time)
#x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)
plt.plot(time,x)
plt.show()
plt.plot(time,x)
plt.xlim([0,2])
plt.show()

f, Pxx_den = signal.welch(x, fs, nperseg=512,scaling='density')
print (f[np.where(Pxx_den==np.max(Pxx_den))])
print (np.max(Pxx_den))
print (simps(Pxx_den,f))

plt.semilogy(f, Pxx_den)
#plt.ylim([0.5e-3, 100])
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [V**2/Hz]')
#plt.show()

f, Pxx_den = signal.welch(x, fs, nperseg=1024,scaling='density')
print (f[np.where(Pxx_den==np.max(Pxx_den))])
print (np.max(Pxx_den))
print (simps(Pxx_den,f))

plt.semilogy(f, Pxx_den)
#plt.ylim([0.5e-3, 100])
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [V**2/Hz]')
plt.show()

# f, Pxx_spec = signal.welch(x, fs, 'flattop', 2048, scaling='spectrum')
# plt.figure()
# plt.semilogy(f, np.sqrt(Pxx_spec))
# plt.xlabel('frequency [Hz]')
# plt.ylabel('Linear spectrum [V RMS]')
# plt.show()

print ('Fiinish!')


