"""
author:yunnaidan
@time: 2018/10/29
@file: pir41event.py
"""
import os
import numpy as np
import pandas as pd
from obspy import UTCDateTime
import multiprocessing
import time
import warnings

from dyntripy.utils import _gen_f_win_list, _f_column_name, load_pi


def _sum_pi(column_name, tar_df):
    interrupt = False
    pi_mean_list = []
    for col in column_name:
        if col not in tar_df.columns:
            raise ValueError('The target frequency range %s not in the database.' % col)
        pi_list = [float(i) for i in tar_df[col].values]
        pi_list = np.array(pi_list)

        # If there is a data gap, we will through this distant event.
        if len(pi_list[pi_list == 0]) != 0:
            interrupt = True

        # Translate sum to mean, so that the time normalization is not needed.
        pi_mean = np.mean(pi_list)
        pi_mean_list.append(pi_mean)

    pi = np.sum(pi_mean_list)

    return interrupt, pi


def _sum_pi_4_f_ranges(sum_f_ranges, f_step, tar_df_b, tar_df_e):
    pi_b_list = []
    pi_e_list = []
    pi_ratio_log_list = []
    interrupt = []
    for ranges in sum_f_ranges:
        f_min = int(ranges[0])
        f_max = int(ranges[1])
        column_name = _f_column_name(f_min, f_step, f_max)
        interrupt_b, pi_b = _sum_pi(column_name, tar_df_b)
        interrupt_e, pi_e = _sum_pi(column_name, tar_df_e)
        pi_ratio = pir(pi_b, pi_e)

        pi_ratio_log = str(np.log10(pi_ratio))
        pi_b = str(pi_b)
        pi_e = str(pi_e)

        pi_b_list.append(pi_b)
        pi_e_list.append(pi_e)
        pi_ratio_log_list.append(pi_ratio_log)
        interrupt.append(interrupt_b)
        interrupt.append(interrupt_e)
    return pi_b_list, pi_e_list, pi_ratio_log_list, interrupt


def pir(pi_b, pi_e):
    # The modify in function _sum_pi makes this time normalization not needed.
    return pi_e / pi_b


def pir41event(
        full_sta,
        event,
        tb_b,
        tb_e,
        te_b,
        te_e,
        f_min,
        f_max,
        f_bin,
        f_step,
        pi_database_folder,
        out_file,
        lock
):
    try:
        file_exist, pi_data_frame = load_pi(
            full_sta, tb_b, te_e, pi_database_folder)

        sum_f_ranges = _gen_f_win_list(f_min, f_bin, f_max)
        sum_f_range_names = _f_column_name(f_min, f_bin, f_max)

        if file_exist:
            pi_data_frame['time'] = [UTCDateTime(t) for t in pi_data_frame['time']]
            tar_df_b = pi_data_frame[(pi_data_frame['time'] >= tb_b) & (pi_data_frame['time'] <= tb_e)]
            tar_df_e = pi_data_frame[(pi_data_frame['time'] >= te_b) & (
                    pi_data_frame['time'] <= te_e)]

            if len(tar_df_b) == 0 or len(tar_df_e) == 0:
                warnings.warn('No pi data of %s for %s!' % (str(event), full_sta))
                pi_b_list = [''] * len(sum_f_ranges)
                pi_e_list = [''] * len(sum_f_ranges)
                pi_ratio_log_list = [''] * len(sum_f_ranges)
            else:
                pi_b_list, pi_e_list, pi_ratio_log_list, interrupt = \
                    _sum_pi_4_f_ranges(sum_f_ranges, f_step, tar_df_b, tar_df_e)

                if True in interrupt:
                    warnings.warn('Interruptions exist around %s for %s!'
                                  % (str(event), full_sta))

        else:
            warnings.warn('No pi data of %s for %s!' % (str(event), full_sta))
            pi_b_list = [''] * len(sum_f_ranges)
            pi_e_list = [''] * len(sum_f_ranges)
            pi_ratio_log_list = [''] * len(sum_f_ranges)

        lock.acquire()
        with open(out_file, 'a+') as f:
            f.write(
                '\n'.join(
                    [
                        ','.join([
                            str(event),
                            f_range_name,
                            pi_b_list[f_range_name_j],
                            pi_e_list[f_range_name_j],
                            pi_ratio_log_list[f_range_name_j]
                        ])
                        for f_range_name_j, f_range_name in enumerate(sum_f_range_names)
                    ]
                )
                + '\n'
            )
        lock.release()

    except Exception as err_msg:
        print(event, full_sta)
        print(err_msg)

    return None


def pir41event_run(
        row,
        f_step,
        PIdatabase_path,
        out_path,
        out_file_prefix,
        lock1,
        lock2
):
    try:
        event = UTCDateTime(row['time'])
        full_sta = row['sta_name']

        tb_b = UTCDateTime(row['Tb_Begin'])
        tb_e = UTCDateTime(row['Tb_End'])
        te_b = UTCDateTime(row['Te_Begin'])
        te_e = UTCDateTime(row['Te_End'])

        f_min = row['f_min']
        f_max = row['f_max']
        f_bin = row['f_bin']
        if len(str(f_min).split('.')) > 1 \
                or len(str(f_max).split('.')) > 1 \
                or len(str(f_bin).split('.')) > 1:
            raise ValueError('The min and max values of frequency must be integer.')
        if f_bin % f_step != 0:
            raise ValueError(
                'The frequency interval of integration must be an integer '
                'multiple of that in the database.'
            )

        pi_database_folder = os.path.join(
            PIdatabase_path,
            full_sta
        )
        out_file = os.path.join(
            out_path,
            out_file_prefix + '_' + full_sta + '.csv'
        )

        lock1.acquire()
        if not os.path.exists(out_file):
            with open(out_file, 'a') as f:
                f.write(','.join(['time', 'f_range', 'Ib', 'Ie', 'ratio']))
                f.write('\n')
        lock1.release()

        with open(out_file, 'r') as f:
            calculated_events = f.readlines()
            calculated_time = [line.split(',')[0] for line in calculated_events[1:]]
            if str(event) not in calculated_time:
                pir41event(
                    full_sta,
                    event,
                    tb_b,
                    tb_e,
                    te_b,
                    te_e,
                    f_min,
                    f_max,
                    f_bin,
                    f_step,
                    pi_database_folder,
                    out_file,
                    lock2,
                )
            else:
                print('PIR of %s for %s exists!' % (str(event), full_sta))
    except Exception as err_msg:
        print(row['time'], row['sta_name'])
        print(err_msg)
    return None


def run_pir_parallel(
        f_step,
        catalog,
        PIdatabase_path,
        out_path,
        out_file_prefix,
        cores
):
    lock1 = multiprocessing.Manager().Lock()
    lock2 = multiprocessing.Manager().Lock()

    event_df = pd.read_csv(catalog)

    pool = multiprocessing.Pool(processes=cores)
    tasks = [(row, f_step, PIdatabase_path, out_path, out_file_prefix, lock1, lock2)
             for row_index, row in event_df.iterrows()]
    # row = event_df.iloc[1509]
    # pir41event_run(row, f_step, PIdatabase_path, out_path, out_file_prefix, lock1, lock2)

    print('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(pir41event_run, tasks, chunksize=1)
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


