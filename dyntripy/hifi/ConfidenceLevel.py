"""
author:yunnaidan
@time: 2018/11/23
@file: ConfidenceLevel.py
"""
import os
import re
import time
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import warnings

from dyntripy.plot_utils import distribution
from dyntripy.utils import _f_column_name


def _remove_abnormal(lst, method='Gaussion'):
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


def _random_cl(
        tele_pir,
        bg_pir_list,
        threshold,
        plot_result=True,
        times=10,
        npertime=30,
        seed=1234,
):
    np.random.seed(seed)
    if len(bg_pir_list) < npertime:
        cl = 'Insufficient RBs'
        mean_cl = ''
        std_cl = ''
        tri_flag = ''
        fig = None
        ax = None
    else:
        cl = []
        fig, ax = plt.subplots()
        for time_i in range(times):
            random_index = np.random.randint(0, len(bg_pir_list), npertime)
            random_bg_list = bg_pir_list[random_index]

            random_bg_list_filtered = _remove_abnormal(random_bg_list, method='Gaussion')
            mean, std = stats.norm.fit(random_bg_list_filtered)
            background_norm = stats.norm(loc=mean, scale=std)
            cl.append(background_norm.cdf(tele_pir))

            if plot_result is not None:
                distribution(
                    random_bg_list_filtered,
                    tele_pir,
                    ax=ax
                )

        mean_cl, std_cl = stats.norm.fit(cl)

        if mean_cl >= threshold:
            tri_flag = str(int(1))
        else:
            tri_flag = str(int(0))

    return cl, mean_cl, std_cl, tri_flag, fig, ax


def cl41event(
        event,
        full_sta,
        sum_f_range_names,
        matched_ratio_df,
        threshold,
        out_file,
        figure_out_folder,
        lock
):
    try:
        cl_list = []
        mean_cl_list = []
        std_cl_list = []
        tri_flag_list = []
        for f_range in sum_f_range_names:
            figure_out_file = os.path.join(
                figure_out_folder,
                str(event) + str(f_range) + '.pdf'
            )

            tar_pir_row = matched_ratio_df[
                (matched_ratio_df['time'] == event) & (matched_ratio_df['f_range'] == f_range)
                ].iloc[0]
            if len(tar_pir_row) == 0:
                warnings.warn('No matched PIR of %s for %s!' % (event, full_sta))
                cl = ''
                mean_cl = ''
                std_cl = ''
                tri_flag = ''
            else:
                tele_pir = tar_pir_row['remote_earthquake_ratio']
                if not np.isnan(tele_pir):
                    tele_pir = float(tele_pir)
                    str_bg_pir_list = re.split(' +', tar_pir_row['background_day_ratio'])
                    bg_pir_list = np.array(str_bg_pir_list, dtype='float64')
                    if figure_out_file is not None:
                        plot_result = True
                    else:
                        plot_result = False
                    cl, mean_cl, std_cl, tri_flag, fig, ax = _random_cl(
                        tele_pir,
                        bg_pir_list,
                        threshold,
                        plot_result=plot_result,
                    )
                    if figure_out_file is not None:
                        ax.set_title(str(event) + ' ' + str(f_range) + '\n'
                                     + 'CL: mean = %.3f' % mean_cl + ' ' + 'std = %.3f' % std_cl)
                        fig.savefig(figure_out_file)
                        plt.close()

                    cl = ' '.join([str(i) for i in cl])
                    mean_cl = str(mean_cl)
                    std_cl = str(std_cl)
                else:
                    cl = 'No RE'
                    mean_cl = ''
                    std_cl = ''
                    tri_flag = ''

            cl_list.append(cl)
            mean_cl_list.append(mean_cl)
            std_cl_list.append(std_cl)
            tri_flag_list.append(tri_flag)

        lock.acquire()
        with open(out_file, 'a') as f:
            f.write(
                '\n'.join(
                    [
                        ','.join([
                            str(event),
                            str(f_range),
                            cl_list[j],
                            str(mean_cl_list[j]),
                            str(std_cl_list[j]),
                            str(tri_flag_list[j])
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


def cl41event_run(
        row,
        matched_ratio_path,
        out_path,
        figure_out_path,
        threshold,
        lock,
):
    try:
        event = UTCDateTime(row['time'])
        full_sta = row['sta_name']

        f_min = row['f_min']
        f_max = row['f_max']
        f_bin = row['f_bin']
        sum_f_range_names = _f_column_name(f_min, f_bin, f_max)

        cal = True
        matched_ratio_file = os.path.join(
            matched_ratio_path,
            'matched_' + full_sta + '.csv')
        if not os.path.exists(matched_ratio_file):
            matched_ratio_df = None
            warnings.warn('No matched PIR file of %s!' % full_sta)
            cal = False
        else:
            matched_ratio_df = pd.read_csv(matched_ratio_file)

        if cal:
            lock.acquire()
            out_file = os.path.join(out_path, 'CL_' + full_sta + '.csv')
            if not os.path.exists(out_file):
                with open(out_file, 'a') as f:
                    f.write(','.join(['time', 'f_range', 'cl', 'cl_mean', 'cl_std', 'triggering']))
                    f.write('\n')
            if figure_out_path != 'None':
                figure_out_folder = os.path.join(figure_out_path, full_sta)
                if not os.path.exists(figure_out_folder):
                    os.makedirs(figure_out_folder)
            else:
                figure_out_folder = None
            lock.release()

            with open(out_file, 'r') as f:
                calculated_events = f.readlines()
                calculated_time = [line.split(',')[0] for line in calculated_events[1:]]
                if str(event) not in calculated_time:
                    cl41event(
                        event,
                        full_sta,
                        sum_f_range_names,
                        matched_ratio_df,
                        threshold,
                        out_file,
                        figure_out_folder,
                        lock
                    )
    except Exception as err_msg:
        print(row['time'], row['sta_name'])
        print(err_msg)
    return None


def run_cl_parallel(
        catalog,
        matched_ratio_path,
        out_path,
        figure_out_path,
        threshold,
        cores
):
    event_df = pd.read_csv(catalog)

    lock = multiprocessing.Manager().Lock()

    pool = multiprocessing.Pool(processes=cores)
    tasks = [(row, matched_ratio_path, out_path, figure_out_path, threshold, lock)
             for row_i, row in event_df.iterrows()]

    # row = event_df.iloc[0]
    # cl41event_run(
    #     row,
    #     matched_ratio_path,
    #     out_path,
    #     figure_out_path,
    #     threshold,
    #     lock
    # )

    print('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(cl41event_run, tasks, chunksize=1)
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
