import numpy as np
from obspy import UTCDateTime
from scipy.stats import norm
from scipy.integrate import simps
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

from dyntrigger.utils.basic_utils import psd


def plot_psd(tele_tr, tb_b, tb_e, te_b, te_e, time_segment, out_file=None, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    def mean_psd(
            time_begin,
            time_end,
            tr_data,
            fs,
            time_segment):
        no = round((time_end - time_begin) / time_segment)
        pxx_all_summary = []
        for n in range(no):
            b_index = round((time_begin + n * time_segment) * fs)
            e_index = round((time_begin + (n + 1) * time_segment) * fs)
            data = tr_data[b_index:e_index]
            pxx_all, f_all = psd(data, fs)
            pxx_all_summary.append(pxx_all)
        pxx_all_summary = np.array(pxx_all_summary)
        mean_pxx_all = np.mean(pxx_all_summary, axis=0)
        return f_all, mean_pxx_all

    tr = tele_tr.copy()
    tr.detrend('linear')
    tr.detrend('constant')
    
    real_tb_b = UTCDateTime(tb_b.year, tb_b.month, tb_b.day) + int(
        (tb_b.hour * 3600 + tb_b.minute * 60 + tb_b.second) / time_segment) * time_segment
    real_tb_e = UTCDateTime(tb_e.year, tb_e.month, tb_e.day) + int(
        (tb_e.hour * 3600 + tb_e.minute * 60 + tb_e.second) / time_segment) * time_segment
    real_te_b = UTCDateTime(te_b.year, te_b.month, te_b.day) + (int(
        (te_b.hour * 3600 + te_b.minute * 60 + te_b.second) / time_segment) + 1) * time_segment
    real_te_e = UTCDateTime(te_e.year, te_e.month, te_e.day) + (int(
        (te_e.hour * 3600 + te_e.minute * 60 + te_e.second) / time_segment) + 1) * time_segment

    fs = tr.stats.sampling_rate
    abs_data_begin = tr.stats.starttime

    tb_b_sec = float('%.2f' % (real_tb_b - abs_data_begin))
    tb_e_sec = float('%.2f' % (real_tb_e - abs_data_begin))
    te_b_sec = float('%.2f' % (real_te_b - abs_data_begin))
    te_e_sec = float('%.2f' % (real_te_e - abs_data_begin))

    f_all_before, mean_pxx_all_before = mean_psd(
        tb_b_sec, tb_e_sec, tr, fs, time_segment)
    f_all_after, mean_pxx_all_after = mean_psd(
        te_b_sec, te_e_sec, tr, fs, time_segment)

    ax.plot(f_all_before, mean_pxx_all_before, '-r')
    ax.plot(f_all_after, mean_pxx_all_after, '-b')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(['$\mathregular{PSD_b}$',
               '$\mathregular{PSD_e}$'],
              frameon=False,
              loc='upper right',
              prop={'style': 'italic'})

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Frequency ' + '$\mathregular{(Hz)}$')
    ax.set_ylabel('Power Spectral Density ' +
                  '$\mathregular{((nm/s)^2/(Hz\cdot s))}$')

    if out_file is not None:
        plt.savefig(out_file)

    return ax


def plot_bg_pir(
        tele_pir,
        background_pirs,
        out_file=None,
        ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    ax.hist(background_pirs, density=True, alpha=0.5, color='gray')

    mean, std = norm.fit(background_pirs)
    bg_norm = norm(loc=mean, scale=std)
    x = np.linspace(bg_norm.ppf(0.001), bg_norm.ppf(0.999), 100)
    y_pdf = bg_norm.pdf(x)
    y_cdf = bg_norm.cdf(x)

    ax.plot(x, y_pdf, '-k')
    # ax.plot(x, y_cdf, '-k')
    ymin, ymax = ax.get_ylim()
    ax.plot([mean, mean], [ymin, ymax - 0.1], '--k')
    ax.plot([tele_pir, tele_pir], [ymin, ymax - 0.1], '--', color='tab:orange')

    ax.set_xlabel('Logarithmic Ratio')
    ax.set_ylabel('Probability Density')

    if out_file is not None:
        plt.savefig(out_file)

    return None

