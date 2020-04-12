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
import logging

logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def remove_abnormal(lst, method='Gaussion'):
    """
    :param lst: Input data of list type.
    :param method: The available methods contains 'IQR', 'Gaussion'
    :return:
    """
    if method == 'IQR':
        Percentile = np.percentile(lst, [0, 25, 50, 75, 100])
        IQR = Percentile[3] - Percentile[1]
        UpLimit = Percentile[3] + IQR * 1.5
        DownLimit = Percentile[1] - IQR * 1.5

        list_filtered = lst[DownLimit <= lst <= UpLimit]

    elif method == 'Gaussion':
        mean_value = np.nanmean(lst)
        var_value = np.sqrt(np.nanvar(lst))
        list_filtered = lst[abs(lst - mean_value) <= 3 * var_value]

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
        background_pir_file,
        pir_column_name,
        day_window,
        out_file):

    tele_df = pd.read_csv(telesesmic_catalog)
    tele_list = [UTCDateTime(i) for i in tele_df['time'].values]

    background_pir_df = pd.read_csv(background_pir_file)

    tele_background_pir_df = pd.DataFrame(
        columns=['time', 'background_PIRs'])

    for tele_index, tele in enumerate(tele_list):
        print('Associate data of %s' % str(tele))
        tele_background_pir_df.loc[tele_index, 'time'] = tele

        tar_background_list = get_tar_background(tele, day_window)
        bg_pir_list = background_pir_df[background_pir_df['time'].isin(
            tar_background_list)][pir_column_name].values

        # Remove the abnormal background event's HighPID data.
        bg_pir_list_filter = remove_abnormal(
            bg_pir_list, method='Gaussion')
        tele_background_pir_df.loc[tele_index,
                                   'background_PIRs'] = bg_pir_list_filter

    tele_background_pir_df.to_csv(out_file, index=False)

    return None
