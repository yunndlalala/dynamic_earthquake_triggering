"""
@version:
author:yunnaidan
@time: 2018/11/27
@file: visual_teleseism.py
@function: Some visualization function about teleseism.
"""
import matplotlib.pyplot as plt
import matplotlib
plt.rc('font',family='Times New Roman')
# matplotlib.rcParams.update({'font.size': 1})
# plt.rc('xtick', labelsize=6)
# plt.rc('ytick', labelsize=6)
from scipy.signal import welch
from scipy.integrate import simps
import os
import numpy as np
import pandas as pd
import obspy
import matplotlib.pyplot as plt
from pylab import *
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
from scipy.signal import spectrogram
from dyntrigger.hifi.ConfidenceLevel import read_associated_background_PIR,load_PIR,view
def psd_interal(pxx_all,f,f_win):
    f_min=f_win[0]
    f_max=f_win[1]
    f_index=np.where((f>= f_min) & (f<=f_max))
    f_tar=f[f_index]
    pxx_tar=pxx_all[f_index]


    PI=simps(pxx_tar,f_tar)
    PI=np.sum(pxx_tar)
    return PI
def spectrogram_plot(event_tr, ax):
    tr = event_tr.copy()
    tr.detrend('linear')
    tr.detrend('constant')
    t1 = float('%.2f' % tr.stats.sac['t1'])
    t2 = float('%.2f' % tr.stats.sac['t2'])
    t3 = float('%.2f' % tr.stats.sac['t3'])
    t0 = float('%.2f' % ((t1 - (t3 - t1)) / 10))
    data_begin = tr.stats.sac['b']
    fs = tr.stats.sampling_rate
    t0_index = int((t0 - data_begin) * fs)
    t3_index = int((t3 - data_begin) * fs)
    cal_tb=t0-60
    cal_te=t3+60
    cal_tb_index=t0_index-6000
    cal_te_index = t3_index +6000
    data = tr.data[cal_tb_index:cal_te_index]
    f_min = 0.1
    f_max = 40
    nfft = 512
    nperseg = nfft
    #noverlap = nfft * 0.75
    #f, t, Sxx = spectrogram(data, fs, nfft=nfft, nperseg=nperseg, noverlap=noverlap)
    f, t, Sxx = spectrogram(data, fs, nfft=nfft, nperseg=nperseg)
    t=t+cal_tb
    img = ax.pcolormesh(t, f, np.log10(abs(Sxx / np.max(Sxx))) * 10, vmin=-160, vmax=0, cmap='jet')
    #img = ax.pcolormesh(t, f, np.log10(abs(Sxx / np.max(Sxx))) * 10, cmap='jet')

    ax.set_ylim(f_min, f_max)
    ax.set_xlim([t0, t3])
    #fig.colorbar(img,ax=ax,orientation='horizontal')

    return img
def waveform_plot(tr, ax,downsample=10):
    tr_copy=tr.copy()
    tr_copy.detrend('linear')
    tr_copy.detrend('constant')
    t1 =float('%.2f' % tr_copy.stats.sac['t1'])
    t2=float('%.2f' % tr_copy.stats.sac['t2'])
    t3 =float('%.2f' % tr_copy.stats.sac['t3'])
    t0 = float('%.2f'%((t1 - (t3 - t1)) / 10))
    data_begin = tr_copy.stats.sac['b']
    fs = tr_copy.stats.sampling_rate
    t0_index = int((t0 - data_begin) * fs)
    t1_index = int((t1 - data_begin) * fs)
    t2_index = int((t2 - data_begin) * fs)
    t3_index = int((t3 - data_begin) * fs)
    data= tr_copy.data[t0_index:t3_index]
    data_downsample=data[np.arange(0,len(data),downsample)]
    ax.plot(np.arange(0,len(data),downsample)*(1/fs)+t0, data_downsample, 'k', linewidth=0.5)
    ax.plot([t1,t1],[np.min(data),np.max(data)],'--r')
    ax.plot([t2, t2], [np.min(data), np.max(data)], '--b')
    #科学计数法
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(xfmt)

    ax.set_xlim([t0, t3])

def waveformFilter_plot(tr, ax,b=None,e=None,downsample=2,freq = 25):
    tr_copy=tr.copy()
    tr_copy.detrend('linear')
    tr_copy.detrend('constant')
    data_begin = float('%.2f'%tr_copy.stats.sac['b'])
    fs = tr_copy.stats.sampling_rate

    tr_copy.filter('highpass', freq=freq)
    if b==None and e==None:
        t1 =float('%.2f' % tr_copy.stats.sac['t1'])
        t2=float('%.2f' % tr_copy.stats.sac['t2'])
        t3 =float('%.2f' % tr_copy.stats.sac['t3'])
        t0 = float('%.2f'%((t1 - (t3 - t1)) / 10))

        t0_index = int((t0 - data_begin) * fs)
        t1_index = int((t1 - data_begin) * fs)
        t2_index = int((t2 - data_begin) * fs)
        t3_index = int((t3 - data_begin) * fs)
    else:
        t0=b
        t3=e
        t0_index = int((t0 - data_begin) * fs)
        t3_index = int((t3 - data_begin) * fs)

    data= tr_copy.data[t0_index:t3_index+1]
    data_downsample = data[np.arange(0, len(data), downsample)]
    ax.plot(np.arange(0, len(data), downsample)*(1/fs)+t0, data_downsample, 'k',linewidth=0.5)
    #ax.plot([t1,t1],[np.min(data),np.max(data)],'--k')
    #ax.plot([t2, t2], [np.min(data), np.max(data)], '-r')

    #ax.text('Highpass > %sHz'%(freq),verticalalignment='top',horizontalalignment='right')

    # text(0.9, 0.9, 'Highpass > %sHz'%(freq),
    #      horizontalalignment='right',
    #      verticalalignment='center',
    #      transform=ax.transAxes)
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
    ax.yaxis.set_major_formatter(xfmt)

    ax.set_xlim([t0,t3])


if __name__ == "__main__":
    teleseism_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/removeRES/GDXB/events/20100404224042.36/NC.GDXB.HHZ'
    event_name="2010-04-04T22:40:42.36Z"
    st=obspy.read(teleseism_file)
    teleseism_tr=st[0]
    sta=teleseism_tr.stats.station
    chn=teleseism_tr.stats.channel
    windowBefore_begin = teleseism_tr.stats.sac['t1']-18000
    windowBefore_end=teleseism_tr.stats.sac['t1']
    windowEvent_begin=teleseism_tr.stats.sac['t2']
    windowEvent_end=teleseism_tr.stats.sac['t3']

    fig=plt.figure(figsize=[12,8])

    #ax1=fig.add_subplot(321)
    ax1 = fig.add_axes([0.05, 0.66, 0.45, 0.25])
    waveform_plot(teleseism_tr,ax1)
    ax1.set_xticklabels([])
    ax1.set_xlabel('')

    ax2 = fig.add_axes([0.05, 0.37, 0.45, 0.25])
    waveformFilter_plot(teleseism_tr,ax2)
    ax2.set_xticklabels([])
    ax2.set_xlabel('')

    ax3 = fig.add_axes([0.05, 0.08, 0.45, 0.25])
    img=spectrogram_plot(teleseism_tr,ax3)
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.8,wspace=0.1, hspace=0.0)
    #plt.subplots_adjust(wspace=0.28, hspace=0.13)
    #ax4=fig.add_subplot(222)
    ax4=fig.add_axes([0.61, 0.54, 0.34, 0.38])
    powerSpectralDensity_plot(teleseism_tr, windowBefore_begin, windowBefore_end, windowEvent_begin, windowEvent_end,
                              ax4)
    #ax5=fig.add_subplot(224)
    ax5 = fig.add_axes([0.61, 0.08, 0.34, 0.38])
    datapath='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Result/Geysers/vel5_vel2_timewindow/fixed_before_window/25_35frequency_range'
    tele_HighPID_file=os.path.join(datapath,'tele_HighPIR_GDXB.csv')
    sortedHighPID_outfile=os.path.join(datapath,'sorted_tele_backgroundHighPIRList_GDXB.txt')
    tele='2010-04-04T22:40:42.36Z'
    backgroundDistribution_plot(tele_HighPID_file, sortedHighPID_outfile, tele, ax5)
    plt.suptitle(event_name + ' Station : %s.%s' % (sta, chn))
    plt.subplots_adjust(wspace=0.28, hspace=0.13)

    ax30 = fig.add_axes([0.51, 0.08, 0.013, 0.24])
    fig.colorbar(img, cax=ax30)

    plt.savefig(
        '/Users/yunnaidan/Project/Research/Dynamic_triggering/Documents/MINE/HighPIR_paper/Figures_Tables/Figure2_20100404/PSD_distribution.pdf',dpi=200)

    fig.show()

    print ('Finish!')

