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
import warnings
import time

from dyntripy.utils import _f_column_name


def _get_tar_background(tele, daywindow):
    tar_bg_list = []
    days = np.arange(-1 * daywindow[0], daywindow[1] + 1, 1)
    for day in days[days != 0]:
        tar_bg = tele + timedelta(days=int(day))
        tar_bg_list.append(str(tar_bg))
    return tar_bg_list


def associate(
        event,
        full_sta,
        sum_f_range_names,
        day_window,
        tele_pir_df,
        background_pir_df,
        pir_column_name,
        out_file,
        lock,
):
    try:
        bg_list = []
        tele_list = []
        for f_range in sum_f_range_names:
            tar_background_list = _get_tar_background(event, day_window)
            tar_bg_pir_df = background_pir_df[
                (background_pir_df['time'].isin(tar_background_list)) &
                (background_pir_df['f_range'] == f_range)]
            tar_bg_pir_list = tar_bg_pir_df[pir_column_name].values
            tar_bg_pir_list = tar_bg_pir_list[~np.isnan(tar_bg_pir_list)]
            tar_bg_pir_list_str = ' '.join([str(i) for i in tar_bg_pir_list])

            tar_tele_pir_df = tele_pir_df[
                       (tele_pir_df['time'] == str(event)) &
                       (tele_pir_df['f_range'] == f_range)]
            tar_tele_pir = tar_tele_pir_df[pir_column_name].values[0]
            if np.isnan(tar_tele_pir):
                tar_tele_pir = ''

            bg_list.append(tar_bg_pir_list_str)
            tele_list.append(tar_tele_pir)

        lock.acquire()
        with open(out_file, 'a') as f:
            f.write(
                '\n'.join(
                    [
                        ','.join([
                            str(event),
                            f_range,
                            str(tele_list[j]),
                            bg_list[j]
                        ])
                        for j, f_range in enumerate(sum_f_range_names)
                    ]
                )
                + '\n'
            )
        lock.release()

    except Exception as err_msg:
        print(event, full_sta)
        print(err_msg)

    return None


def associate_run(
        row,
        day_window,
        pir_column_name,
        re_path,
        rb_path,
        out_path,
        lock1,
        lock2
):
    try:
        event = UTCDateTime(row['time'])
        full_sta = row['sta_name']

        f_min = row['f_min']
        f_max = row['f_max']
        f_bin = row['f_bin']
        sum_f_range_names = _f_column_name(f_min, f_bin, f_max)

        cal = True
        tele_pir_file = os.path.join(
            re_path, 'RE_' + full_sta + '.csv')
        if not os.path.exists(tele_pir_file):
            warnings.warn('No PIR file of %s!' % full_sta)
            cal = False
        background_pir_file = os.path.join(
            rb_path, 'RB_' + full_sta + '.csv')
        if not os.path.exists(background_pir_file):
            warnings.warn('No background PIR file of %s!' % full_sta)
            cal = False

        if cal:
            tele_pir_df = pd.read_csv(tele_pir_file)
            background_pir_df = pd.read_csv(background_pir_file)

            out_file = os.path.join(
                out_path, 'matched_' + full_sta + '.csv')

            lock1.acquire()
            if not os.path.exists(out_file):
                with open(out_file, 'a') as f:
                    f.write(','.join(['time', 'f_range', 'remote_earthquake_ratio', 'background_day_ratio']))
                    f.write('\n')
            lock1.release()

            with open(out_file, 'r') as f:
                calculated_events = f.readlines()
                calculated_time = [line.split(',')[0] for line in calculated_events[1:]]
                if str(event) not in calculated_time:
                    associate(
                        event,
                        full_sta,
                        sum_f_range_names,
                        day_window,
                        tele_pir_df,
                        background_pir_df,
                        pir_column_name,
                        out_file,
                        lock2
                    )
                else:
                    print('Matched PIR of %s for %s exists!' % (str(event), full_sta))
    except Exception as err_msg:
        print(row['time'], row['sta_name'])
        print(err_msg)

    return None


def run_associate_parallel(
        telesesmic_catalog,
        pir_column_name,
        day_window,
        re_path,
        rb_path,
        out_path,
        cores
):
    lock1 = multiprocessing.Manager().Lock()
    lock2 = multiprocessing.Manager().Lock()

    catalog_df = pd.read_csv(telesesmic_catalog)

    pool = multiprocessing.Pool(processes=cores)
    tasks = [(row, day_window, pir_column_name, re_path, rb_path, out_path, lock1, lock2)
             for row_i, row in catalog_df.iterrows()]
    # row = catalog_df.iloc[0]
    # associate_run(row, day_window, pir_column_name, re_path, rb_path, out_path, lock1, lock2)

    print('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(associate_run, tasks, chunksize=1)
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
    print('\n')
