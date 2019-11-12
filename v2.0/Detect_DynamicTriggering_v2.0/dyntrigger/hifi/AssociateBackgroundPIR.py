"""
@version:
author:yunnaidan
@time: 2018/09/12
@file: AssociateBackgroundPIR.py
@function:
"""

import numpy as np
import pandas as pd
import os
from obspy.core import UTCDateTime
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

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

def get_uniqe(list):
    out_list = [list[0]]
    out_list_index = [0]
    for i in range(1, len(list)):
        if list[i] not in out_list:
            out_list.append(list[i])
            out_list_index.append(i)
    return out_list


def load_PIR(PIR_file, column_name):
    dataframe = pd.read_csv(PIR_file)
    events = get_uniqe(dataframe['event'].values)
    event_PIR_dict = {}
    for i in range(len(events)):
        event = events[i]
        PIR= dataframe[dataframe['event']== event][column_name].values[0]
        event_PIR_dict[event] = PIR
    return event_PIR_dict


def get_tar_background(tele, daywindow):
    tar_bg_list = []
    days = np.arange(-1 * (daywindow/2), (daywindow/2) + 1, 1)
    for day in days[days != 0]:
        tar_bg = tele + timedelta(days=int(day))
        tar_bg_list.append(str(tar_bg))
    return tar_bg_list


def save_data(tele_background_PIR_dict, outfile):
    with open(outfile, 'w') as f:
        for key in tele_background_PIR_dict.keys():
            f.write(str(key))
            for i in tele_background_PIR_dict[key]:
                f.write(',' + str(i))
            f.write('\n')


def PIR_associate(
        telesesmic_catalog,
        background_PIR_file,
        PIR_column_name,
        dayWindow,
        out_file):

    tele_dataframe = pd.read_csv(telesesmic_catalog)
    tele_list = [UTCDateTime(i) for i in tele_dataframe['time'].values]

    background_PIR_dict = load_PIR(
        background_PIR_file, PIR_column_name)

    tele_background_PIR_dict = {}
    for tele in tele_list:
        print('Associate data of teleseism : %s' % str(tele))
        tar_background_list = get_tar_background(tele, dayWindow)

        bgPIR_list = []
        for tar_bg in tar_background_list:
            # It is possible that no data for one background event.
            if tar_bg in background_PIR_dict.keys():
                # bgHighPID_list.append(backgroundTele_HighPID_dict[tar_bgTele])
                # 先进行log操作再去掉异常值，实验表明如果反过来结果相差会比较大。因此在记录的文件中也只能记录log之后的
                bgPIR_list.append(background_PIR_dict[tar_bg])
        # Remove the abnormal background event's HighPID data.
        _, bgPIR_list_filter = remove_abnormal(
            bgPIR_list, method='Gaussion')
        tele_background_PIR_dict[str(tele)] = bgPIR_list_filter
    save_data(tele_background_PIR_dict, out_file)
    return tele_background_PIR_dict


def plot_sort_data(
        teleCatalog_file,
        backgroundTele_HighPID_file,
        daywindow,
        figure_outpath):

    dataframe = pd.read_csv(teleCatalog_file, encoding="gb2312")
    tele_list = [UTCDateTime(i) for i in dataframe['time'].values]

    dataframe_background = pd.read_csv(
        backgroundTele_HighPID_file, encoding="gb2312")
    background_tele_list = [
        UTCDateTime(i) for i in dataframe_background['event'].values]

    meanPIR_winRatio_list=[]

    for tele in tele_list:
        fig = plt.figure()
        # tele=UTCDateTime('2009-08-05T09:13:12.30Z')
        print('Sort data of teleseism : %s' % str(tele))
        tar_backgroundTele_list = get_tar_background(tele, daywindow)

        PIR_list = []
        PI_list = []

        for tar_bgTele in tar_backgroundTele_list:
            # It is possible that no data for one background event.
            if UTCDateTime(tar_bgTele) in background_tele_list:
                str_PI_b_list = dataframe_background[dataframe_background['event']
                                                     == tar_bgTele]['PI_b_list'].values
                str_PI_e_list = dataframe_background[dataframe_background['event']
                                                     == tar_bgTele]['PI_e_list'].values
                PI_b_list = [float(i) for i in str_PI_b_list[0].split(';')]
                PI_e_list = [float(i) for i in str_PI_e_list[0].split(';')]
                PIR = np.log10(np.mean(PI_e_list) / np.mean(PI_b_list))

                PIR_list.append(PIR)
                PI_list.append([PI_b_list, PI_e_list])
        filterd_index, filterd_data = remove_abnormal(
            PIR_list, method='Gaussion')

        PI_b_mean_list = []
        PI_e_mean_list = []
        for index in filterd_index:
            PI_b_mean_list.append(np.log10(np.mean(PI_list[index][0])))
            PI_e_mean_list.append(np.log10(np.mean(PI_list[index][1])))
            plt.semilogy(PI_list[index][0], '-k')
            plt.semilogy(range(len(PI_list[index][0]), len(
                PI_list[index][0]) + len(PI_list[index][1])), PI_list[index][1], '-',color='gray')

        for i in range(len(PI_b_mean_list)):
            x_b_median=len(PI_list[0][0]) / 2
            x_e_median=len(PI_list[0][0]) +(len(PI_list[0][1]) /2)
            plt.plot([x_b_median,x_e_median],[10**PI_b_mean_list[i],10**PI_e_mean_list[i]],'o')
        xlim_min=0
        xlim_max=len(PI_list[0][0]) + len(PI_list[0][1]) + 1
        plt.plot([xlim_min,xlim_max],
                 [10**np.mean(PI_b_mean_list),10**np.mean(PI_b_mean_list)],'--b')
        plt.plot([xlim_min,xlim_max],
                 [10**np.mean(PI_e_mean_list),10**np.mean(PI_e_mean_list)],'--r')

        # median_type = dict(color='darkgreen')
        # mean_type = dict(color='black')
        # plt.boxplot([PI_b_mean_list,PI_e_mean_list],
        #             positions=[150,301],
        #             widths=10,
        #             meanline=True,
        #             showmeans=True,
        #             whis=3,
        #             medianprops=median_type,
        #             meanprops=mean_type)
        plt.xlim([xlim_min,xlim_max])
        plt.title(tele)
        plt.savefig(os.path.join(figure_outpath, str(tele) + '.png'))
    return None


def test_windowRatio(
        teleCatalog_file,
        backgroundTele_HighPID_file,
        daywindow):
    dataframe = pd.read_csv(teleCatalog_file, encoding="gb2312")
    tele_list = [UTCDateTime(i) for i in dataframe['time'].values]

    dataframe_background = pd.read_csv(
        backgroundTele_HighPID_file, encoding="gb2312")
    background_tele_list = [
        UTCDateTime(i) for i in dataframe_background['event'].values]

    meanPIR_winRatio_list = []

    for tele in tele_list:
        # tele=UTCDateTime('2009-08-05T09:13:12.30Z')
        print('Sort data of teleseism : %s' % str(tele))
        tar_backgroundTele_list = get_tar_background(tele, daywindow)

        PIR_list = []
        PI_list = []

        for tar_bgTele in tar_backgroundTele_list:
            # It is possible that no data for one background event.
            if UTCDateTime(tar_bgTele) in background_tele_list:
                str_PI_b_list = dataframe_background[dataframe_background['event']
                                                     == tar_bgTele]['PI_b_list'].values
                str_PI_e_list = dataframe_background[dataframe_background['event']
                                                     == tar_bgTele]['PI_e_list'].values
                PI_b_list = [float(i) for i in str_PI_b_list[0].split(';')]
                PI_e_list = [float(i) for i in str_PI_e_list[0].split(';')]
                PIR = np.log10(np.mean(PI_e_list) / np.mean(PI_b_list))

                PIR_list.append(PIR)
                PI_list.append([PI_b_list, PI_e_list])
        filterd_index, filterd_data = remove_abnormal(
            PIR_list, method='Gaussion')
        meanPIR=np.mean(filterd_data)
        winRatio=np.log10(len(PI_list[0][1])/len(PI_list[0][0]))
        meanPIR_winRatio_list.append([meanPIR,winRatio])
    for i in range(len(meanPIR_winRatio_list)):
        plt.plot(2 * i, meanPIR_winRatio_list[i][0], 'o', color='black')
        plt.plot(2 * i, meanPIR_winRatio_list[i][1], 'o', color='red')
        plt.plot([2 * i, 2 * i], [meanPIR_winRatio_list[i][0], meanPIR_winRatio_list[i][1]], '-', color='gray')
        plt.plot(2 * i, meanPIR_winRatio_list[i][0] - meanPIR_winRatio_list[i][1], '*', color='black')
    plt.plot([0, 150], [0, 0], '--k')
    plt.show()



if __name__ == '__main__':

    print('finish')
