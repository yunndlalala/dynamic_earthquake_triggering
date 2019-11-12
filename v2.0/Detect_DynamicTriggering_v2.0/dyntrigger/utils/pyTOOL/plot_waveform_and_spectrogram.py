"""
@version:
author:yunnaidan
@time: 2018/10/01
@file: plot_waveform_and_spectrogram.py
@function:
"""
import os
import pandas as pd
import numpy as np
import obspy
from obspy.signal.tf_misfit import cwt
from obspy.imaging.cm import obspy_sequential
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
from scipy.signal import spectrogram
def plot_wave_and_spec4oneEvent(event_path,sta,chn,time_before=-3600,time_after=3600):
    """
    :param event_path: should be absolute path
    :param sta: ex. 'GDXB'
    :param chn: ex. 'HHZ'
    :param time_before: should be minus
    :param time_after: 
    :return: 
    """
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2])
    ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60], sharex=ax1)
    ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])

    st=obspy.read(os.path.join(event_path,'*'+sta+'*'+chn+'*'))
    tr=st[0]
    tr.detrend('linear')
    tr.detrend('constant')
    #tr.taper(0.5, type='hann')

    t1 =float('%.2f'%tr.stats.sac['t1'])
    t2=float('%.2f'%tr.stats.sac['t2'])
    t3 =float('%.2f'%tr.stats.sac['t3'])
    t0 = float('%.2f'%((t1 - (t3 - t1)) / 10))
    data_begin = tr.stats.sac['b']
    fs = tr.stats.sampling_rate
    t0_index = int((t0 - data_begin) * fs)
    t1_index = int((t1 - data_begin) * fs)
    t2_index = int((t2 - data_begin) * fs)
    t3_index = int((t3 - data_begin) * fs)
    data=tr.data[t0_index:t3_index]
    ax1.plot(np.arange(t0,t3,1/fs), data, 'k')
    ax1.set_xlim([t0,t3])
    f_min = 0.1
    f_max = 40
    #scalogram = cwt(data, 1/fs, 8, f_min, f_max) #计算速度太慢了！！
    nfft = 1024;
    nperseg = nfft;
    noverlap = nfft * 0.75
    f, t, Sxx = spectrogram(data, fs,nfft=nfft,nperseg=nperseg,noverlap=noverlap)
    #cmap = plt.get_cmap('PiYG')
    #norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    img=ax2.pcolormesh(t, f, np.log10(abs(Sxx/np.max(Sxx)))*10,vmin=-150,vmax=0,cmap='jet')

    #spectrogram(x, fs=1.0, window=(‘tukey’, 0.25), nperseg=None, noverlap=None, nfft=None,
    #                         detrend=’constant’, return_onesided = True, scaling =’density’, axis = -1, mode =’psd’)
    # x, y = np.meshgrid(
    #     np.arange(t0, t3, 1 / fs),
    #     np.logspace(np.log10(f_min), np.log10(f_max), scalogram.shape[0]))
    # x, y = np.meshgrid(
    #     np.arange(t0, t3, 1 / fs),
    #     np.arange(f_min, f_max, (f_max-f_min)/scalogram.shape[0]))
    # img=ax2.pcolormesh(x, y, np.log10(abs(scalogram/np.max(scalogram)))*(10), cmap=obspy_sequential)
    ax2.set_xlabel("Time after %s [s]" % tr.stats.starttime)
    ax2.set_ylabel("Frequency [Hz]")
    #ax2.set_yscale('log')
    ax2.set_ylim(f_min, f_max)
    ax2.set_xlim([t0,t3])
    fig.colorbar(img, cax=ax3)


    # amarker=tr.stats.sac['a']
    # sampling_rate=tr.stats.sampling_rate
    # b_time=amarker+time_before
    # e_time=amarker+time_after
    # b_point=b_time*sampling_rate
    # e_point=e_time*sampling_rate
    fig.show()

        
if __name__=="__main__":
    event_path='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/removeRES/GDXB/events/20100404224042.36'
    sta='NC.GDXB'
    chn='HHZ'
    plot_wave_and_spec4oneEvent(event_path, sta, chn, time_before=-3600, time_after=3600)
    #get_catalog(origin_data_file, out_file, before_time_in_s=0, cal_Parrivetime=True)
    print ('Finish!')