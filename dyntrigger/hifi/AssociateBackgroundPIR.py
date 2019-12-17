"""
@version:
author:yunnaidan
@time: 2018/09/12
@file: AssociateBackgroundPIR.py
@function:
"""

import numpy as np
import pandas as pd
from obspy.core import UTCDateTime
from datetime import timedelta


def remove_abnormal(list, method='Gaussion'):
    """
    :param list: Inputdata of list type.
    :param method: The available methods contains 'IQR', 'Gaussion'
    :return:
    """
    if method == 'IQR':
        Percentile = np.percentile(list, [0, 25, 50, 75, 100])
        IQR = Percentile[3] - Percentile[1]
        UpLimit = Percentile[3] + IQR * 1.5
        DownLimit = Percentile[1] - IQR * 1.5

        list_filtered = list[DownLimit <= list <= UpLimit]

    elif method == 'Gaussion':
        mean_value = np.nanmean(list)
        var_value = np.sqrt(np.nanvar(list))
        list_filtered = list[abs(list - mean_value) <= 3 * var_value]

    else:
        raise ValueError('The param of method is not supported!')

    return list_filtered


def get_tar_background(tele, daywindow):
    tar_bg_list = []
    days = np.arange(-1 * (daywindow / 2), (daywindow / 2) + 1, 1)
    for day in days[days != 0]:
        tar_bg = tele + timedelta(days=int(day))
        tar_bg_list.append(str(tar_bg))
    return tar_bg_list


def run_associate(
        telesesmic_catalog,
        background_PIR_file,
        PIR_column_name,
        dayWindow,
        out_file):

    tele_df = pd.read_csv(telesesmic_catalog)
    tele_list = [UTCDateTime(i) for i in tele_df['time'].values]

    background_PIR_df = pd.read_csv(background_PIR_file)

    tele_background_PIR_df = pd.DataFrame(
        columns=['time', 'background_PIRs'])

    for tele_index, tele in enumerate(tele_list):
        # print('Associate data of teleseism : %s' % str(tele))
        tele_background_PIR_df.loc[tele_index, 'time'] = tele

        tar_background_list = get_tar_background(tele, dayWindow)
        bgPIR_list = background_PIR_df[background_PIR_df['event'].isin(
            tar_background_list)][PIR_column_name].values

        # Remove the abnormal background event's HighPID data.
        bgPIR_list_filter = remove_abnormal(
            bgPIR_list, method='Gaussion')
        tele_background_PIR_df.loc[tele_index,
                                   'background_PIRs'] = bgPIR_list_filter

    tele_background_PIR_df.to_csv(out_file, index=False)

    return tele_background_PIR_df
