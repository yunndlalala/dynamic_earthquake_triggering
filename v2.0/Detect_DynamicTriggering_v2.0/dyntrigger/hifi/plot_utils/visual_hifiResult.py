"""
@version:
author:yunnaidan
@time:
@file: visual_hifiResult.py
@function: Some visualization function for the results of the HiFi method.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman', size=8)
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
import matplotlib.cbook as cbook
from obspy import UTCDateTime
from dyntrigger.hifi.ConfidenceLevel import read_associated_background_PIR, load_PIR, view
from scipy.signal import welch
from scipy.integrate import simps
import obspy
from dyntrigger.hifi.PIdatabase import Gf


def get_uniqe(list):
    out_list = [list[0]]
    out_list_index = [0]
    for i in range(1, len(list)):
        if list[i] not in out_list:
            out_list.append(list[i])
            out_list_index.append(i)
    return out_list, out_list_index


def remove_abnormal(list, method='Gaussion'):
    """
    :param list: Inputdata of list type.
    :param method: The available methods contains 'Boxplot', 'Gaussion'
    :return:
    """
    if method == 'Boxplot':
        stats = cbook.boxplot_stats(list, whis=10)
        outers = stats[0]['fliers']
        list_filter = []
        index_filter = []
        for i in range(len(list)):
            item = list[i]
            if item not in outers:
                list_filter.append(item)
                index_filter.append(i)
    elif method == 'Gaussion':
        mean_value = np.mean(list)
        var_value = np.var(list) ** (0.5)
        list_filter = []
        index_filter = []
        for i in range(len(list)):
            item = list[i]
            if abs(item - mean_value) <= 3 * var_value:
                list_filter.append(item)
                index_filter.append(i)
    else:
        raise ValueError('The param of method is not supported!')
    return index_filter, list_filter


'''
    The follow functions are used to displace the results of multiple stations.
'''


def HighPID4oneteleseism_plot(
        teleseism_name,
        HighPID_dataframe,
        reabnormal=False,
        show_object='all'):
    """
    :param teleseism_name:
    :param HighPID_dataframe:
    :param reabnormal:
    :param show_object: three choices contains 'all'/'filter'/'raw'
    :return:
    """
    def run_plot(sta_list, HighPID_list, ax, title_name='None'):
        ax.scatter([i for i in range(len(HighPID_list))], HighPID_list, s=120)
        ax.set_xticks([i for i in range(len(sta_list))])
        ax.set_xticklabels(sta_list)
        ax.set_title(title_name)
        plt.xlabel('Stations')
        plt.ylabel('Power Integral Difference/(nm/s)^2')
        return ax
    [teleseism_columnName, sta_columnName,
        HighPID_columnName] = HighPID_dataframe.columns.values.tolist()
    tar_dataframe = HighPID_dataframe[HighPID_dataframe[teleseism_columnName]
                                      == teleseism_name]
    tar_HighPID = tar_dataframe[HighPID_columnName].values
    tar_sta = tar_dataframe[sta_columnName].values
    if reabnormal:
        index_filter, tar_HighPID_filter = remove_abnormal(
            tar_HighPID, method='Gaussion')
        tar_sta_filter = []
        for i in range(len(tar_dataframe)):
            if i in index_filter:
                tar_sta_filter.append(tar_dataframe.iloc[i][sta_columnName])
    fig = plt.figure(figsize=[10, 8])
    if show_object == 'all':
        ax1 = fig.add_subplot(211)
        run_plot(tar_sta, tar_HighPID, ax1, name='original')
        ax2 = fig.add_subplot(212)
        run_plot(
            tar_sta_filter,
            tar_HighPID_filter,
            ax2,
            name='after remove the outliers')
    elif show_object == 'raw':
        ax = fig.add_subplot(111)
        run_plot(tar_sta, tar_HighPID, ax, name='original')
    elif show_object == 'filter':
        ax1 = fig.add_subplot(111)
        run_plot(
            tar_sta_filter,
            tar_HighPID_filter,
            ax1,
            name='after remove the outliers')
    else:
        raise ValueError('The input parameter "show" is wrong!')

    return fig


def HighPID4oneStation_plot(
        station_name,
        HighPID_dataframe,
        show_object='all'):
    """
    :param station_name:
    :param HighPID_dataframe:
    :param show_object: contains three kind 'all'/'withnoise'/'nonoise'
    :return:
    """

    if show_object == 'all':
        [teleseism_columnName, sta_columnName, HighPID_nonoise_columnName,
            HighPID_withnoise_columnName] = HighPID_dataframe.columns.values.tolist()
        tar_dataframe = HighPID_dataframe[HighPID_dataframe[sta_columnName]
                                          == station_name]
        tar_withnoise_HighPID = tar_dataframe[HighPID_withnoise_columnName].values
        tar_nonoise_HighPID = tar_dataframe[HighPID_nonoise_columnName].values
    elif show_object == 'withnoise':
        [teleseism_columnName, sta_columnName,
            HighPID_withnoise_columnName] = HighPID_dataframe.columns.values.tolist()
        tar_dataframe = HighPID_dataframe[HighPID_dataframe[sta_columnName]
                                          == station_name]
        tar_withnoise_HighPID = tar_dataframe[HighPID_withnoise_columnName].values
    elif show_object == 'nonoise':
        [teleseism_columnName, sta_columnName,
            HighPID_nonoise_columnName] = HighPID_dataframe.columns.values.tolist()
        tar_dataframe = HighPID_dataframe[HighPID_dataframe[sta_columnName]
                                          == station_name]
        tar_nonoise_HighPID = tar_dataframe[HighPID_nonoise_columnName].values
    else:
        raise ValueError('The input parameter "show" is wrong!')
    tar_teleseism_datetime = [
        UTCDateTime(i) for i in tar_dataframe[teleseism_columnName].values]
    tar_teleseism = ['/'.join([str(tele.year).zfill(2),
                               str(tele.month).zfill(2),
                               str(tele.day).zfill(2),
                               str(tele.hour).zfill(2)]) for tele in tar_teleseism_datetime]

    fig = plt.figure(figsize=[18, 10])
    ax = fig.add_subplot(111)
    if show_object == 'all':
        ax.scatter(
            tar_teleseism_datetime,
            tar_withnoise_HighPID,
            c='',
            edgecolors='black',
            marker='o',
            s=400)
        ax.scatter(
            tar_teleseism_datetime,
            tar_nonoise_HighPID,
            marker='+',
            c='black',
            s=400)
    elif show_object == 'withnoise':
        ax.scatter(
            tar_teleseism_datetime,
            tar_withnoise_HighPID,
            c='',
            edgecolors='black',
            marker='o',
            s=400)
    elif show_object == 'nonoise':
        ax.scatter(
            tar_teleseism_datetime,
            tar_nonoise_HighPID,
            marker='+',
            c='black',
            s=400)

    ax.set_xticks([i for i in tar_teleseism_datetime])
    ax.set_xticklabels(tar_teleseism)
    plt.xticks(rotation=90, fontsize=18)
    plt.yticks(fontsize=20)
    plt.title('The HighPID of GDXB station for every teleseism', fontsize=20)
    plt.xlabel('Teleseism', fontsize=20)
    plt.ylabel('Power Integral Difference/(nm/s)^2', fontsize=20)

    return fig


'''
    The follow functions are used to displace the results of single stations.
'''


def background_distribution_plot(
        tele_PIR_file,
        associated_file,
        tele,
        ax, **kwags):

    tele_background_PIR_dict = read_associated_background_PIR(
        associated_file)
    tele_PIR_dict = load_PIR(tele_PIR_file, 'PIRatio_log')

    tele = str(UTCDateTime(tele))
    PIR_list = tele_background_PIR_dict[tele]
    tele_PIR = float(tele_PIR_dict[tele])

    samples = view(PIR_list, tele_PIR, ax, **kwags)

    return samples


def confidenceLevel_plot(datafile, ax):
    dataframe = pd.read_csv(datafile)
    for i in range(len(dataframe)):
        event_datetime = UTCDateTime(dataframe.iloc[i]['tele'])
        PIRatio = float(dataframe.iloc[i]['confidence_level_value'])
        trigger_refer = int(dataframe.iloc[i]['reference_trigger'])
        beta_value = float(dataframe.iloc[i]['beta_value'])

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
    ax.set_xticks([UTCDateTime(2000, 1, 1),
                   UTCDateTime(2001, 1, 1),
                   UTCDateTime(2002, 1, 1),
                   UTCDateTime(2003, 1, 1),
                   UTCDateTime(2004, 1, 1),
                   UTCDateTime(2005, 1, 1),
                   UTCDateTime(2006, 1, 1),
                   UTCDateTime(2007, 1, 1),
                   UTCDateTime(2008, 1, 1),
                   UTCDateTime(2009, 1, 1),
                   UTCDateTime(2010, 1, 1),
                   UTCDateTime(2011, 1, 1),
                   UTCDateTime(2012, 1, 1),
                   UTCDateTime(2013, 1, 1)])
    ax.set_xticklabels(['2000', '2001',
                        '2002', '2003',
                        '2004', '2005',
                        '2006', '2007',
                        '2008', '2009',
                        '2010', '2011',
                        '2012', '2013'])
    ax.plot([UTCDateTime(2000, 1, 1), UTCDateTime(
        2013, 1, 1)], [0.025, 0.025], '--r')
    ax.text(UTCDateTime(2000, 1, 1), 0.1, 'Threshold: 0.025', size=14)
    ax.set_xlabel('Date')
    ax.set_ylabel('Significance  Level')
    return None


def confidenceLevel_vs_betaValue_plot(datafile, ax):
    dataframe = pd.read_csv(datafile)
    for i in range(len(dataframe)):
        event_datetime = UTCDateTime(dataframe.iloc[i]['tele'])
        confidence_level = float(dataframe.iloc[i]['confidence_level_mean'])
        trigger_refer = int(dataframe.iloc[i]['reference_trigger'])
        beta_value = float(dataframe.iloc[i]['beta_value'])
        error = float(dataframe.iloc[i]['confidence_level_var'])

        #ax.errorbar(confidence_level, beta_value,xerr=error, fmt='o', markersize=5,ecolor='gray', elinewidth=1, capsize=2)
        ax.errorbar(
            confidence_level,
            beta_value,
            xerr=error,
            fmt='o',
            markersize=5,
            ecolor='gray',
            color='white',
            markeredgecolor='k',
            elinewidth=2,
            capsize=3)

    # ax.plot([0, 1], [2, 2], '--r')
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

    def mean_psd(time_begin, time_end, tr_data, fs, Gf_parameters,time_segment):  # time_begin and time_end (s)
        N = int((time_end - time_begin) / time_segment)
        pxx_all_summary = []
        for n in range(N):
            data = tr_data[int((time_begin + n * time_segment) * fs)
                               :int((time_begin + (n + 1) * time_segment) * fs)]
            pxx_all, f_all = psd(data, fs)
            sensivity, normalizing, zeros, poles = Gf_parameters
            Gf_list = [Gf(sensivity, normalizing, zeros, poles, f)
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
        Tb_B_time, Tb_E_time, teleseism_tr, fs, Gf_parameters,time_segment)
    f_all_after, mean_pxx_all_after = mean_psd(
        Te_B_time, Te_E_time, teleseism_tr, fs, Gf_parameters,time_segment)


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


if __name__ == '__main__':
    tele_origin = '2009-08-05T09:13:12.30Z'
    catalog = '/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Geysers/Data/Geysers/vel5_vel2_timewindow/fixed_before_window/25_35frequency_range/EQ_info_GeysersGDXB.csv'
    data_file = '/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Geysers/Data/Raw_data/20090805/NC.GDXB.HHZ'
    st = obspy.read(data_file)
    teleseism_tr = st[0]
    # A_time=UTCDateTime('2010-04-04T22:42:55.420781Z')
    # B_time = UTCDateTime('2010-04-04T22:44:07.140000Z')
    # E_time = UTCDateTime('2010-04-04T22:49:14.310000Z')
    A_time = UTCDateTime('2009-08-05T09:16:04.916744Z')
    B_time = UTCDateTime('2009-08-05T09:17:41.220000Z')
    E_time = UTCDateTime('2009-08-05T09:24:24.600000Z')

    sensivity = 1.006000e+09
    zeros = [0,
             0,
             complex(-4.618140e+02, -4.290780e+02),
             complex(-4.618140e+02, 4.290780e+02),
             complex(-1.946970e+02, 0),
             complex(-1.514870e+01, 0)]
    poles = [complex(-3.681800e-02, -3.692300e-02),
             complex(-3.681800e-02, 3.692300e-02),
             complex(-1.023970e+04, -2.725020e+03),
             complex(-1.023970e+04, 2.725020e+03),
             complex(-9.512740e+03, -1.146990e+04),
             complex(-9.512740e+03, 1.146990e+04),
             complex(-4.545250e+02, 0.000000e+00),
             complex(-3.981840e+02, 0.000000e+00),
             complex(-9.294710e+01, -3.889920e+02),
             complex(-9.294710e+01, 3.889920e+02),
             complex(-1.545030e+01, 0.000000e+00)]
    normalizing = 9.482280e+18
    Gf_parameters = [sensivity, normalizing, zeros, poles]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    tick_width = 0.5
    tick_length = 2
    powerSpectralDensity_plot(
        teleseism_tr,
        A_time,
        B_time,
        E_time,
        Gf_parameters,
        ax)

    ax.plot([25, 25], [10 ** (-10), 10 ** (5)], '--k')
    ax.plot([35, 35], [10 ** (-10), 10 ** (5)], '--k')
    ax.tick_params(
        which='both',
        direction='in',
        length=tick_length,
        width=tick_width)
    ax.set_xlim([1, 40])
    #ax.set_ylim([10 ** (-10), 10 ** (10)])
    ax.set_xlabel(
        'Frequency ' +
        '$\mathregular{(Hz)}$')
    fig.show()
    print('Finish')
