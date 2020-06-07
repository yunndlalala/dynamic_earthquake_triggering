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
import logging


logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def load_pi(full_sta, tb_b, te_e, pi_database_path):
    b_day = ''.join([str(tb_b.year),
                     str(tb_b.month).zfill(2),
                     str(tb_b.day).zfill(2)])
    e_day = ''.join([str(te_e.year),
                     str(te_e.month).zfill(2),
                     str(te_e.day).zfill(2)])
    file_exist = False
    if b_day == e_day:
        tar_file = os.path.join(
            pi_database_path,
            full_sta +
            '_' +
            str(b_day) +
            '.csv')
        if os.path.exists(tar_file):
            file_exist = True
            pi_data_frame = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    full_sta +
                    '_' +
                    str(b_day) +
                    '.csv'))
        else:
            pi_data_frame = None
    else:
        tar_file1 = os.path.join(
            pi_database_path, full_sta + '_' + str(b_day) + '.csv')
        tar_file2 = os.path.join(
            pi_database_path, full_sta + '_' + str(e_day) + '.csv')
        if os.path.exists(tar_file1) and os.path.exists(tar_file2):
            file_exist = True
            pi_data_frame1 = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    full_sta +
                    '_' +
                    str(b_day) +
                    '.csv'))
            pi_data_frame2 = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    full_sta +
                    '_' +
                    str(e_day) +
                    '.csv'))
            pi_data_frame = pi_data_frame1.append(
                pi_data_frame2, ignore_index=True)
        else:
            pi_data_frame = None
    return file_exist, pi_data_frame


def f_column_name(f_min, f_max, step=5):
    bin_no = int((f_max - f_min) / step)
    column_name = []
    for i in range(bin_no):
        name = str(f_min + i * step) + '-' + str(f_min + (i + 1) * step)
        column_name.append(name)
    return column_name


def sum_pi(column_name, tar_df):
    interrupt = False
    pi_mean_list = []
    for col in column_name:
        pi_list = [float(i) for i in tar_df[col].values]
        pi_list = np.array(pi_list)
        # If there is a data gap larger than 60s, we will through this distant
        # event.
        if len(pi_list[pi_list == 0]) != 0:
            interrupt = True
            break
        else:
            # Translate sum to mean, so that the time normalization is not needed.
            # This modify has been applied and tested, the changes of results
            # are really small.
            pi_mean = np.mean(pi_list)
            pi_mean_list.append(pi_mean)
    if interrupt:
        pi = 0
    else:
        pi = np.sum(pi_mean_list)

    list4plot = np.zeros(len(tar_df))
    for col in column_name:
        pi_list = [float(i) for i in tar_df[col].values]
        list4plot += pi_list

    return interrupt, pi, list4plot


def pir(pi_b, pi_e):
    # The modify in function sum_pi makes this time normalization not needed.
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
        f_step,
        pi_database_path,
        out_file):

    file_exist, pi_data_frame = load_pi(
        full_sta, tb_b, te_e, pi_database_path)

    if file_exist:
        pi_data_frame['time'] = [UTCDateTime(t)
                                for t in pi_data_frame['time']]
        tar_df_b0 = pi_data_frame[(pi_data_frame['time'] >= tb_b)
                                 & (pi_data_frame['time'] <= tb_e)]
        tar_df_b = pi_data_frame.iloc[tar_df_b0.index - 1]
        tar_df_e = pi_data_frame[(pi_data_frame['time'] >= te_b) & (
                pi_data_frame['time'] <= te_e)]

        column_name = f_column_name(f_min, f_max, step=f_step)
        interrupt_b, pi_b, list4plot_b = sum_pi(column_name, tar_df_b)
        interrupt_e, pi_e, list4plot_e = sum_pi(column_name, tar_df_e)

        if interrupt_b or interrupt_e:
            interrupt_final = True
        else:
            interrupt_final = False
            pi_ratio = pir(pi_b, pi_e)
            pi_ratio_log = np.log10(pi_ratio)

        if interrupt_final:
            logging.warning(
                'There is interrupt in this event %s!' %
                str(event))
            pi_b = None
            pi_e = None
            pi_ratio_log = None
    else:
        logging.warning('There is no data for this event %s!' % str(event))
        pi_b = None
        pi_e = None
        pi_ratio_log = None

    with open(out_file, 'a') as f:
        f.write(','.join([str(event)[:-4] + 'Z', str(pi_b), str(
            pi_e), str(pi_ratio_log)]))
        f.write('\n')

    return None


def run_pir_parallel(f_step, catalog, PIdatabase_path, full_sta, out_file, cores):
    if os.path.exists(out_file):
        os.remove(out_file)
    with open(out_file, 'a') as f:
        f.write(','.join(['time', 'Ib', 'Ie', 'ratio']))
        f.write('\n')

    event_df = pd.read_csv(catalog)
    pool = multiprocessing.Pool(processes=cores)
    tasks = []
    for row_index, row in event_df.iterrows():
        event = UTCDateTime(event_df.iloc[row_index]['time'])
        print('Prepare ratio task %s' % event, end='\r')

        tb_b = UTCDateTime(event_df.iloc[row_index]['Tb_Begin'])
        tb_e = UTCDateTime(event_df.iloc[row_index]['Tb_End'])
        te_b = UTCDateTime(event_df.iloc[row_index]['Te_Begin'])
        te_e = UTCDateTime(event_df.iloc[row_index]['Te_End'])

        f_min = event_df.iloc[row_index]['f_min']
        f_max = event_df.iloc[row_index]['f_max']

        tasks.append((full_sta,
                      event,
                      tb_b,
                      tb_e,
                      te_b,
                      te_e,
                      f_min,
                      f_max,
                      f_step,
                      PIdatabase_path,
                      out_file))

        # pir41event(
        #     full_sta,
        #     event,
        #     tb_b,
        #     tb_e,
        #     te_b,
        #     te_e,
        #     f_min,
        #     f_max,
        #     f_step,
        #     PIdatabase_path,
        #     out_file)
    print ('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(pir41event, tasks, chunksize=1)
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


