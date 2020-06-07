"""
@version:
author:yunnaidan
@time: 2018/09/12
@file: AssociateBackgroundPIR.py
@function:
"""
import os
import numpy as np
import pandas as pd
from obspy.core import UTCDateTime
from datetime import timedelta
import multiprocessing
import time
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
    days = np.arange(-1 * daywindow[0], daywindow[1] + 1, 1)
    for day in days[days != 0]:
        tar_bg = tele + timedelta(days=int(day))
        tar_bg_list.append(str(tar_bg)[:-4] + 'Z')
    return tar_bg_list


def associate(tele, day_window, tele_pir_df, background_pir_df, pir_column_name, out_file):
    tar_background_list = get_tar_background(tele, day_window)
    bg_pir_list_str = background_pir_df[background_pir_df['time'].isin(
        tar_background_list)][pir_column_name].values
    bg_pir_list_str = bg_pir_list_str[bg_pir_list_str != 'None']
    bg_pir_list = bg_pir_list_str.astype(np.float64)
    # Remove the abnormal data.
    bg_pir_list_filter = remove_abnormal(
        bg_pir_list, method='Gaussion')
    bg_pir_list_filter_str = ' '.join([str(i) for i in bg_pir_list_filter])

    tele_pir = str(tele_pir_df[tele_pir_df['time'] == str(tele)[:-4] + 'Z'][pir_column_name].values[0])

    with open(out_file, 'a') as f:
        f.write(','.join([str(tele)[:-4] + 'Z', tele_pir, bg_pir_list_filter_str]))
        f.write('\n')

    return None


def run_associate_parallel(
        telesesmic_catalog,
        tele_pir_file,
        background_pir_file,
        pir_column_name,
        day_window,
        out_file,
        cores):
    if os.path.exists(out_file):
        os.remove(out_file)
    with open(out_file, 'a') as f:
        f.write(','.join(['time', 'remote_earthquake_ratio', 'background_day_ratio']))
        f.write('\n')

    tele_df = pd.read_csv(telesesmic_catalog)
    tele_list = [UTCDateTime(i) for i in tele_df['time'].values]

    background_pir_df = pd.read_csv(background_pir_file)
    tele_pir_df = pd.read_csv(tele_pir_file)

    pool = multiprocessing.Pool(processes=cores)
    tasks = []
    for tele_index, tele in enumerate(tele_list):
        print('Prepare matching task %s' % str(tele), end='\r')
        # associate(tele, day_window, tele_pir_df, background_pir_df, pir_column_name, out_file)
        tasks.append((tele, day_window, tele_pir_df, background_pir_df, pir_column_name, out_file))
    print ('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(associate, tasks, chunksize=1)
    # close() & join() is necessary
    pool.close()
    # simple progress bar
    while (True):
        remaining = rs._number_left
        print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
              end='\r')  # '\r' means remove the last line
        if (rs.ready()):
            break
        time.sleep(0.5)

    pool.join()
    print ('\n')
