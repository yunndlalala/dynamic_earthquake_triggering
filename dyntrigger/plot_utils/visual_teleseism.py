"""
author:yunnaidan
@time: 2018/11/27
@file: visual_teleseism.py
"""
from pylab import *
from scipy.signal import spectrogram


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
    cal_tb = t0 - 60
    cal_te = t3 + 60
    cal_tb_index = t0_index - 6000
    cal_te_index = t3_index + 6000
    data = tr.data[cal_tb_index:cal_te_index]
    f_min = 0.1
    f_max = 40
    nfft = 512
    nperseg = nfft
    #noverlap = nfft * 0.75
    #f, t, Sxx = spectrogram(data, fs, nfft=nfft, nperseg=nperseg, noverlap=noverlap)
    f, t, Sxx = spectrogram(data, fs, nfft=nfft, nperseg=nperseg)
    t = t + cal_tb
    img = ax.pcolormesh(t, f, np.log10(abs(Sxx / np.max(Sxx)))
                        * 10, vmin=-160, vmax=0, cmap='jet')
    #img = ax.pcolormesh(t, f, np.log10(abs(Sxx / np.max(Sxx))) * 10, cmap='jet')

    ax.set_ylim(f_min, f_max)
    ax.set_xlim([t0, t3])
    # fig.colorbar(img,ax=ax,orientation='horizontal')

    return img


def waveform_plot(tr, ax, downsample=10):
    tr_copy = tr.copy()
    tr_copy.detrend('linear')
    tr_copy.detrend('constant')
    t1 = float('%.2f' % tr_copy.stats.sac['t1'])
    t2 = float('%.2f' % tr_copy.stats.sac['t2'])
    t3 = float('%.2f' % tr_copy.stats.sac['t3'])
    t0 = float('%.2f' % ((t1 - (t3 - t1)) / 10))
    data_begin = tr_copy.stats.sac['b']
    fs = tr_copy.stats.sampling_rate
    t0_index = int((t0 - data_begin) * fs)
    t1_index = int((t1 - data_begin) * fs)
    t2_index = int((t2 - data_begin) * fs)
    t3_index = int((t3 - data_begin) * fs)
    data = tr_copy.data[t0_index:t3_index]
    data_downsample = data[np.arange(0, len(data), downsample)]
    ax.plot(np.arange(0, len(data), downsample) * (1 / fs) +
            t0, data_downsample, 'k', linewidth=0.5)
    ax.plot([t1, t1], [np.min(data), np.max(data)], '--r')
    ax.plot([t2, t2], [np.min(data), np.max(data)], '--b')
    # 科学计数法
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(xfmt)

    ax.set_xlim([t0, t3])


def waveformFilter_plot(tr, ax, b=None, e=None, downsample=2, freq=25):
    tr_copy = tr.copy()
    tr_copy.detrend('linear')
    tr_copy.detrend('constant')
    data_begin = float('%.2f' % tr_copy.stats.sac['b'])
    fs = tr_copy.stats.sampling_rate

    tr_copy.filter('highpass', freq=freq)
    if b is None and e is None:
        t1 = float('%.2f' % tr_copy.stats.sac['t1'])
        t2 = float('%.2f' % tr_copy.stats.sac['t2'])
        t3 = float('%.2f' % tr_copy.stats.sac['t3'])
        t0 = float('%.2f' % ((t1 - (t3 - t1)) / 10))

        t0_index = int((t0 - data_begin) * fs)
        t1_index = int((t1 - data_begin) * fs)
        t2_index = int((t2 - data_begin) * fs)
        t3_index = int((t3 - data_begin) * fs)
    else:
        t0 = b
        t3 = e
        t0_index = int((t0 - data_begin) * fs)
        t3_index = int((t3 - data_begin) * fs)

    data = tr_copy.data[t0_index:t3_index + 1]
    data_downsample = data[np.arange(0, len(data), downsample)]
    ax.plot(np.arange(0, len(data), downsample) * (1 / fs) +
            t0, data_downsample, 'k', linewidth=0.5)
    # ax.plot([t1,t1],[np.min(data),np.max(data)],'--k')
    #ax.plot([t2, t2], [np.min(data), np.max(data)], '-r')

    #ax.text('Highpass > %sHz'%(freq),verticalalignment='top',horizontalalignment='right')

    # text(0.9, 0.9, 'Highpass > %sHz'%(freq),
    #      horizontalalignment='right',
    #      verticalalignment='center',
    #      transform=ax.transAxes)
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
    ax.yaxis.set_major_formatter(xfmt)

    ax.set_xlim([t0, t3])
