"""
@version:
author:yunnaidan
@time:
@file: visual_hifiResult.py
@function: Some visualization function for the results of the HiFi method.
"""

import pandas as pd
import numpy as np
from obspy import UTCDateTime
from scipy.signal import welch
from scipy.stats import norm

from dyntrigger.hifi.PIdatabase import gf


def bg_pir_plot(
        tele_pir,
        background_pirs,
        ax):
    ax.hist(background_pirs, density=True, alpha=0.5, color='gray')

    mean, std = norm.fit(background_pirs)
    bg_norm = norm(loc=mean, scale=std)
    x = np.linspace(bg_norm.ppf(0.001), bg_norm.ppf(0.999), 100)
    y_pdf = bg_norm.pdf(x)
    y_cdf = bg_norm.cdf(x)

    ax.plot(x, y_pdf, '-k')
    ax.plot(x, y_cdf, '-k')
    ymin, ymax = ax.get_ylim()
    ax.plot([mean, mean], [ymin, ymax - 0.1], '--k')
    ax.plot([tele_pir, tele_pir], [ymin, ymax - 0.1], '--', color='tab:orange')

    return None


def cl_plot(datafile, ax):
    dataframe = pd.read_csv(datafile)
    for i in range(len(dataframe)):
        event_datetime = UTCDateTime(dataframe.iloc[i]['tele'])
        PIRatio = float(dataframe.iloc[i]['confidence_level_value'])
        trigger_refer = int(dataframe.iloc[i]['reference_trigger'])

        if trigger_refer == 0 and PIRatio > 0.025:
            ax.scatter(
                event_datetime,
                PIRatio,
                c='tab:gray',
                edgecolors='tab:gray',
                marker='o',
                s=150)
        if trigger_refer == 0 and PIRatio <= 0.025:
            ax.scatter(
                event_datetime,
                PIRatio,
                c='palegreen',
                edgecolors='black',
                marker='o',
                s=150)
        if trigger_refer == 1 and PIRatio > 0.025:
            ax.scatter(
                event_datetime,
                PIRatio,
                c='lightcoral',
                edgecolors='black',
                marker='o',
                s=150)
        if trigger_refer == 1 and PIRatio <= 0.025:
            ax.scatter(
                event_datetime,
                PIRatio,
                c='olive',
                edgecolors='black',
                marker='o',
                s=150)

    return None


def powerSpectralDensity_plot(
        teleseism_tr,
        Tb_B,
        Tb_E,
        Te_B,
        Te_E,
        Gf_parameters,
        time_segment,
        ax):

    def psd(data, fs):
        f, pxx_all = welch(data,
                           fs,
                           window='hanning',
                           nperseg=512,
                           noverlap=256,
                           nfft=512,
                           detrend=None,
                           return_onesided=True,
                           scaling='density',
                           axis=-1)
        return pxx_all, f

    # time_begin and time_end (s)
    def mean_psd(
            time_begin,
            time_end,
            tr_data,
            fs,
            Gf_parameters,
            time_segment):
        N = int((time_end - time_begin) / time_segment)
        pxx_all_summary = []
        for n in range(N):
            data = tr_data[int((time_begin + n * time_segment) * fs)
                               :int((time_begin + (n + 1) * time_segment) * fs)]
            pxx_all, f_all = psd(data, fs)
            sensivity, normalizing, zeros, poles = Gf_parameters
            Gf_list = [gf(sensivity, normalizing, zeros, poles, f)
                       for f in f_all]
            counts = len(Gf_list)
            pxx_all_removeGf = [(pxx_all[i] / (Gf_list[i] ** 2))
                                * (10 ** 18) for i in range(counts)]
            pxx_all_summary.append(pxx_all_removeGf)
        pxx_all_summary = np.array(pxx_all_summary)
        mean_pxx_all = np.mean(pxx_all_summary, axis=0)
        return f_all, mean_pxx_all

    teleseism_tr.detrend('linear')
    teleseism_tr.detrend('constant')

    fs = teleseism_tr.stats.sampling_rate
    abs_data_begin = teleseism_tr.stats.starttime
    tar_Tb_B = UTCDateTime(Tb_B.year, Tb_B.month, Tb_B.day) + int(
        (Tb_B.hour * 3600 + Tb_B.minute * 60 + Tb_B.second) / time_segment) * time_segment
    tar_Tb_E = UTCDateTime(Tb_E.year, Tb_E.month, Tb_E.day) + int(
        (Tb_E.hour * 3600 + Tb_E.minute * 60 + Tb_E.second) / time_segment) * time_segment
    tar_Te_B = UTCDateTime(Te_B.year, Te_B.month, Te_B.day) + (int(
        (Te_B.hour * 3600 + Te_B.minute * 60 + Te_B.second) / time_segment) + 1) * time_segment
    tar_Te_E = UTCDateTime(Te_E.year, Te_E.month, Te_E.day) + (int(
        (Te_E.hour * 3600 + Te_E.minute * 60 + Te_E.second) / time_segment) + 1) * time_segment

    Tb_B_time = tar_Tb_B - abs_data_begin
    Tb_E_time = tar_Tb_E - abs_data_begin
    Te_B_time = tar_Te_B - abs_data_begin
    Te_E_time = tar_Te_E - abs_data_begin

    f_all_before, mean_pxx_all_before = mean_psd(
        Tb_B_time, Tb_E_time, teleseism_tr, fs, Gf_parameters, time_segment)
    f_all_after, mean_pxx_all_after = mean_psd(
        Te_B_time, Te_E_time, teleseism_tr, fs, Gf_parameters, time_segment)

    ax.plot(f_all_before, mean_pxx_all_before, '-r')
    ax.plot(f_all_after, mean_pxx_all_after, '-b')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(['$\mathregular{PSD_b}$',
               '$\mathregular{PSD_b}$'],
              frameon=False,
              loc='upper left',
              prop={'size': 7,
                    'family': 'Times New Roman',
                    'style': 'italic'})

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    return None