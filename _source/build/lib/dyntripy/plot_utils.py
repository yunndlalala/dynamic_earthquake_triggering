#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2020/08/09
@file: plot_utils.py
"""
import os
import obspy as ob
import numpy as np
from scipy import stats
from scipy.signal import spectrogram
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import warnings

from dyntripy.utils import \
    get_act_windows, select_event_waveforms, load_gf, high_frequency_energy, gf, _gen_f_win_list, _f_column_name


def distribution(
        background_pir_list,
        tele_pir,
        out_file=None,
        ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    mean, std = stats.norm.fit(background_pir_list)
    background_norm = stats.norm(loc=mean, scale=std)
    confidence_level_value = background_norm.cdf(tele_pir)
    x = np.linspace(
        background_norm.ppf(0.001),
        background_norm.ppf(0.999),
        100)
    y_pdf = background_norm.pdf(x)

    ax.hist(
        background_pir_list,
        density=True,
        alpha=0.2,
        # color='gray'
    )
    ax.plot(x, y_pdf, '-', color='gray')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # ax.plot([mean, mean], [ymin, ymax], '-k')
    ax.plot([tele_pir, tele_pir], [ymin, ymax], '--', color='k')
    # ax.text(xmin + 0.8 * (xmax - xmin),
    #         ymin + 0.9 * (ymax - ymin),
    #         'CL = %.3f' % confidence_level_value)
    # ax.set_ylim([ymin, ymax])
    ax.set_xlabel('Logarithmic Ratio')
    ax.set_ylabel('Probability Density')

    if out_file is not None:
        ax_top = ax.twiny()
        ax_top.set_xticks([mean - 3 * std,
                           mean - 2 * std,
                           mean - std,
                           mean,
                           mean + std,
                           mean + 2 * std,
                           mean + 3 * std])
        ax_top.set_xticklabels(
            ['0.001', '0.023', '0.159', '0.500', '0.841', '0.977', '0.999'])
        ax_top.set_xlabel('Confidence Level')
        xlim = ax.get_xlim()
        ax_top.set_xlim(xlim)

        plt.savefig(out_file)
        plt.close()

    return ax


def psds_comparison(
        row,
        waveform_path,
        gf_info_file,
        nperseg,
        noverlap,
        nfft,
        f_step,
        output_path,
):
    try:
        event = row['time']
        full_sta = row['sta_name']
        tb_b = ob.UTCDateTime(row['Tb_Begin'])
        tb_e = ob.UTCDateTime(row['Tb_End'])
        te_b = ob.UTCDateTime(row['Te_Begin'])
        te_e = ob.UTCDateTime(row['Te_End'])
        freq_min = row['f_min']
        freq_max = row['f_max']
        freq_bin = row['f_bin']

        event_st, event_sac_file = select_event_waveforms(full_sta, tb_b, te_e, waveform_path, return_sac_file=True)

        ground_tb_pxx_summary = []
        ground_te_pxx_summary = []
        if event_st is not None:
            for tr_i in range(len(event_st)):
                tr = event_st[tr_i]
                tr.detrend('linear')
                tr.detrend('constant')
                fs = tr.stats.sampling_rate
                starttime = tr.stats.starttime

                sac_file = event_sac_file[tr_i].split('/')[-1]
                _, sensitivity, normalizing, zeros, poles = load_gf(sac_file, gf_info_file)

                f, t, pxx_all = spectrogram(
                    tr.data, fs,
                    window='hanning', nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=None,
                    return_onesided=True, scaling='density', axis=-1, mode='psd'
                )

                t = np.full(len(t), starttime) + t
                tb_index = np.where((t >= tb_b) & (t <= tb_e))[0]
                te_index = np.where((t >= te_b) & (t <= te_e))[0]
                tb_pxx_all = pxx_all[:, tb_index]
                te_pxx_all = pxx_all[:, te_index]

                gf_list = np.array([gf(sensitivity, normalizing, zeros, poles, f_i) for f_i in f])
                ground_tb_pxx = np.array([
                    (tb_pxx_all[i, :] / (gf_list[i] ** 2)) * (10 ** 18)
                    for i in range(len(gf_list))
                ])
                ground_te_pxx = np.array([
                    (te_pxx_all[i, :] / (gf_list[i] ** 2)) * (10 ** 18)
                    for i in range(len(gf_list))
                ])
                if len(ground_tb_pxx_summary) == 0:
                    ground_tb_pxx_summary = ground_tb_pxx
                else:
                    ground_tb_pxx_summary = np.column_stack(ground_tb_pxx_summary, ground_tb_pxx)
                if len(ground_te_pxx_summary) == 0:
                    ground_te_pxx_summary = ground_te_pxx
                else:
                    ground_te_pxx_summary = np.column_stack(ground_te_pxx_summary, ground_te_pxx)
            ground_tb_mean_pxx = np.mean(ground_tb_pxx_summary, axis=1)
            ground_te_mean_pxx = np.mean(ground_te_pxx_summary, axis=1)

            fig, ax = plt.subplots()
            ax.plot(f, ground_tb_mean_pxx, label='$T_b$')
            ax.plot(f, ground_te_mean_pxx, label='$T_e$')
            ax.legend()
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('PSD ((nm/s)^2/Hz)')
            ax.set_title(str(event) + ' ' + full_sta)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            sum_f_ranges = _gen_f_win_list(freq_min, freq_bin, freq_max)
            pi_b_list = []
            pi_e_list = []
            for f_range in sum_f_ranges:
                sub_f_min, sub_f_max = f_range
                pi_b, pi_e = high_frequency_energy(
                    f_b=f,
                    ground_pxx_mean_b=ground_tb_mean_pxx,
                    f_e=f,
                    ground_pxx_mean_e=ground_te_mean_pxx,
                    f_min=sub_f_min,
                    f_max=sub_f_max,
                    f_step=f_step,
                )
                pi_b_list.append(pi_b)
                pi_e_list.append(pi_e)
            ax.text(xmin, ymin,
                    '\n'.join(
                        [
                            str(sum_f_ranges[i][0]) + '-' + str(sum_f_ranges[i][1]) + ' Hz ' +
                            '$I_b$=%i' % int(pi_b_list[i]) + ' ' +
                            '$I_e$=%i' % int(pi_e_list[i]) + ' ' +
                            'ratio=%.2f' % np.log10(pi_e_list[i] / pi_b_list[i])
                            for i in range(len(pi_b_list))]
                    )
                    )
            ax.set_title(str(event) + ' ' + full_sta)
            output_file = os.path.join(
                output_path,
                str(event) + '_' + full_sta + '.png')
            plt.savefig(output_file)
            plt.close()
        else:
            warnings.warn('No pi data of %s.' % full_sta)
    except Exception as err_msg:
        print(event, full_sta)
        print('Error: ', err_msg)
    return None


def spec_wave_highpass(
        row,
        pi_database_path,
        waveform_path,
        sperseg,
        nperseg,
        noverlap,
        nfft,
        wave=True,
        highpass=True,
        output_path=None,
        downsample=1,
        amplitude_unit='Count',
        left=0.13,
        bottom=0.1,
        width=0.24,
        length=0.7,
        bin=0.06,
        before_show_l=0.2,
        after_show_l=0.05,
        cmap='coolwarm',
        vmin=-150,
        vmax=0,
        ax=None,
):
    try:
        ot = ob.UTCDateTime(row['time'])
        full_sta = row['sta_name']
        tb_b = ob.UTCDateTime(row['Tb_Begin'])
        tb_e = ob.UTCDateTime(row['Tb_End'])
        te_b = ob.UTCDateTime(row['Te_Begin'])
        te_e = ob.UTCDateTime(row['Te_End'])
        freq_min = row['f_min']
        freq_max = row['f_max']

        pi_database_folder = os.path.join(pi_database_path, full_sta)
        act_tb_b, act_tb_e, act_te_b, act_te_e = get_act_windows(
                event=ot,
                full_sta=full_sta,
                pi_database_folder=pi_database_folder,
                tb_b=tb_b,
                tb_e=tb_e,
                te_b=te_b,
                te_e=te_e,
                sperseg=sperseg,
        )

        event_st = select_event_waveforms(full_sta, tb_b, te_e, waveform_path)

        if act_tb_b is not None and event_st is not None:
            show_e = act_te_e + (act_te_e - act_te_b) * after_show_l
            show_b = act_tb_e - (act_te_e - act_te_b) * before_show_l

            if ax is None:
                ax1 = plt.axes([left, bottom + 2 * width + 2 * bin, length, width])
                colorbar_ax = plt.axes([left + length + 0.01,
                                        bottom + 2 * width + 2 * bin,
                                        0.03, width])
                ax2 = plt.axes([left, bottom + width + bin, length, width])
                ax3 = plt.axes([left, bottom, length, width])
            else:
                ax1 = ax
                colorbar_ax = None
                ax2 = None
                ax3 = None

            for tr in event_st:
                tr.detrend('linear')
                tr.detrend('constant')
                fs = tr.stats.sampling_rate
                starttime = tr.stats.starttime

                f, t, pxx_all = spectrogram(
                    tr.data, fs,
                    window='hanning', nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=None,
                    return_onesided=True, scaling='density', axis=-1, mode='psd'
                )
                t = np.full(len(t), starttime) + t
                tar_index = np.where((t >= show_b) & (t <= show_e))[0]
                tar_t = t[tar_index]
                tar_pxx_all = pxx_all[:, tar_index]

                plot_t = tar_t - ot - sperseg / 2
                plot_t = np.append(plot_t, plot_t[-1] + sperseg / 2)
                plot_f = f - (f[1] - f[0]) / 2
                plot_f = np.append(plot_f, plot_f[-1] + (f[1] - f[0]))
                img = ax1.pcolormesh(
                    plot_t, plot_f, 10 * np.log10(tar_pxx_all / np.max(tar_pxx_all)), vmin=vmin, vmax=vmax, cmap=cmap
                )
                ax1.plot([show_b - ot, show_e - ot], [freq_min, freq_min], '--k')
                ax1.plot([show_b - ot, show_e - ot], [freq_max, freq_max], '--k')

                if wave:
                    show_tr = tr.slice(show_b, show_e)
                    data = show_tr.data
                    downsample_index = np.arange(0, len(data), downsample)
                    data_downsample = data[downsample_index]
                    time = downsample_index * (1 / fs) + (show_b - ot)
                    ax2.plot(time, data_downsample, 'gray', linewidth=0.5)
                    ax2.plot([act_tb_b - ot, act_tb_b - ot], [np.min(data), np.max(data)], '--k')
                    ax2.plot([act_tb_e - ot, act_tb_e - ot], [np.min(data), np.max(data)], '--k')
                    ax2.plot([act_te_b - ot, act_te_b - ot], [np.min(data), np.max(data)], '-k')
                    ax2.plot([act_te_e - ot, act_te_e - ot], [np.min(data), np.max(data)], '-k')
                if highpass:
                    show_tr = tr.slice(show_b, show_e)
                    show_tr.taper(0.05)
                    show_tr.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
                    data = show_tr.data
                    downsample_index = np.arange(0, len(data), downsample)
                    data_downsample = data[downsample_index]
                    time = downsample_index * (1 / fs) + (show_b - ot)
                    ax3.plot(time, data_downsample, 'k', linewidth=0.5)

            ax1.set_xlim([show_b - ot, show_e - ot])
            ax1.set_ylabel('Frequency (Hz)')

            if wave:
                ax1.set_xticklabels([])
                ax2.set_xlim([show_b - ot, show_e - ot])
                ax2.set_ylabel(amplitude_unit)
                # Scientific notation
                xfmt = ScalarFormatter(useMathText=True)
                xfmt.set_powerlimits((0, 0))
                ax2.yaxis.set_major_formatter(xfmt)
            if highpass:
                ax2.set_xticklabels([])
                ax3.set_xlim([show_b - ot, show_e - ot])
                ax3.set_ylabel(amplitude_unit)
                ax3.set_xlabel('Time (s)')
                xfmt = ScalarFormatter(useMathText=True)
                xfmt.set_powerlimits((0, 0))
                ax3.yaxis.set_major_formatter(xfmt)

            if colorbar_ax is not None:
                cbar = plt.colorbar(img, cax=colorbar_ax)
                cbar.ax.set_xlabel('dB')

            if output_path is not None:
                ax1.set_title(str(ot) + ' ' + full_sta)
                output_file = os.path.join(
                    output_path,
                    str(ot) + '_' + full_sta + '.png')
                plt.savefig(output_file)
                plt.close()
    except Exception as err_msg:
        print(ot, full_sta)
        print('Error: ', err_msg)
    return ax1

