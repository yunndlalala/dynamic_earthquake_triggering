"""
author:yunnaidan
@time: 2018/10/29
@file: pir41event.py
"""
import os
import numpy as np
import pandas as pd
from obspy.core.utcdatetime import UTCDateTime
import logging


logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def load_pi(sta, chn, Tb_B, Te_E, PIdatabase_path):
    B_day = ''.join([str(Tb_B.year),
                     str(Tb_B.month).zfill(2),
                     str(Tb_B.day).zfill(2)])
    E_day = ''.join([str(Te_E.year),
                     str(Te_E.month).zfill(2),
                     str(Te_E.day).zfill(2)])
    file_exist = False
    if B_day == E_day:
        tar_file = os.path.join(
            PIdatabase_path,
            str(sta) + '.' + chn +
            '_' +
            str(B_day) +
            '.csv')
        if os.path.exists(tar_file):
            file_exist = True
            PIdata_frame = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) + '.' + chn +
                    '_' +
                    str(B_day) +
                    '.csv'))
        else:
            PIdata_frame = None
    else:
        tar_file1 = os.path.join(
            PIdatabase_path, str(sta) + '.' + chn + '_' + str(B_day) + '.csv')
        tar_file2 = os.path.join(
            PIdatabase_path, str(sta) + '.' + chn + '_' + str(E_day) + '.csv')
        if os.path.exists(tar_file1) and os.path.exists(tar_file2):
            file_exist = True
            PIdata_frame1 = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) + '.' + chn +
                    '_' +
                    str(B_day) +
                    '.csv'))
            PIdata_frame2 = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) + '.' + chn +
                    '_' +
                    str(E_day) +
                    '.csv'))
            PIdata_frame = PIdata_frame1.append(
                PIdata_frame2, ignore_index=True)
        else:
            PIdata_frame = None
    return file_exist, PIdata_frame


def f_column_name(f_min, f_max, step = 5):
    binNo = int((f_max - f_min) / step)
    column_name = []
    for i in range(binNo):
        name = str(f_min + i * step) + '-' + str(f_min + (i + 1) * step)
        column_name.append(name)
    return column_name


def sum_pi(column_name, tar_df):
    interrupt = False
    PI_mean_list = []
    for col in column_name:
        PI_list = [float(i) for i in tar_df[col].values]
        PI_list = np.array(PI_list)
        # If there is a data gap larger than 60s, we will through this distant
        # event.
        if len(PI_list[PI_list == 0]) != 0:
            interrupt = True
            break
        else:
            # Translate sum to mean, so that the time normalization is not needed.
            # This modify has been applied and tested, the changes of results
            # are really small.
            PI_mean = np.mean(PI_list)
            PI_mean_list.append(PI_mean)
    if interrupt:
        PI = 0
    else:
        PI = np.sum(PI_mean_list)

    list4plot = np.zeros(len(tar_df))
    for col in column_name:
        PI_list = [float(i) for i in tar_df[col].values]
        list4plot += PI_list

    return interrupt, PI, list4plot


def pir(PI_b, PI_e):
    # The modify in function sum_pi makes this time normalization not needed.
    PI_e_2 = PI_e
    PI_b_2 = PI_b
    return PI_e_2 / PI_b_2


def pir41event(
        sta,
        chn,
        event,
        Tb_B,
        Tb_E,
        Te_B,
        Te_E,
        f_min,
        f_max,
        f_step,
        PIdatabase_path,
        out_file):

    file_exist, PIdata_frame = load_pi(
        sta, chn, Tb_B, Te_E, PIdatabase_path)

    if file_exist:
        PIdata_frame['time'] = [UTCDateTime(t)
                                for t in PIdata_frame['time']]
        tar_df_b0 = PIdata_frame[(PIdata_frame['time'] >= Tb_B)
                                 & (PIdata_frame['time'] <= Tb_E)]
        tar_df_b = PIdata_frame.iloc[tar_df_b0.index - 1]
        tar_df_e = PIdata_frame[(PIdata_frame['time'] >= Te_B) & (
            PIdata_frame['time'] <= Te_E)]

        column_name = f_column_name(f_min, f_max, step=f_step)
        interrupt_b, PI_b, list4plot_b = sum_pi(column_name, tar_df_b)
        interrupt_e, PI_e, list4plot_e = sum_pi(column_name, tar_df_e)

        if interrupt_b or interrupt_e:
            interrupt_final = True
            PI_ratio = 0
        else:
            interrupt_final = False
            PI_ratio = pir(PI_b, PI_e)
            PI_ratio_log = np.log10(PI_ratio)

        if interrupt_final:
            logging.warning(
                'There is interrupt in this event %s!' %
                str(event))
        else:
            with open(out_file, 'a') as f:
                f.write(','.join([str(event), str(PI_b), str(
                    PI_e), str(PI_ratio), str(PI_ratio_log)]))
                f.write('\n')

    # If there is no data for the day from Tb to Te, we will through this
    # distant event.
    else:
        logging.warning('There is no data for this event %s!' % str(event))

    return None


def run_pir(catalog, PIdatabase_path, sta, chn, out_file):
    with open(out_file, 'a') as f:
        f.write(','.join(['time', 'PI_b', 'PI_e', 'PIR', 'PIR_log']))
        f.write('\n')

    event_df = pd.read_csv(catalog)
    for row_index, row in event_df.iterrows():
        event = UTCDateTime(event_df.iloc[row_index]['time'])
        print('Calculate pir of %s' % event)

        Tb_B = UTCDateTime(event_df.iloc[row_index]['Tb_B_time'])
        Tb_E = UTCDateTime(event_df.iloc[row_index]['Tb_E_time'])
        Te_B = UTCDateTime(event_df.iloc[row_index]['Te_B_time'])
        Te_E = UTCDateTime(event_df.iloc[row_index]['Te_E_time'])

        f_min = event_df.iloc[row_index]['f_min']
        f_max = event_df.iloc[row_index]['f_max']
        f_step = event_df.iloc[row_index]['f_step']
        pir41event(
            sta,
            chn,
            event,
            Tb_B,
            Tb_E,
            Te_B,
            Te_E,
            f_min,
            f_max,
            f_step,
            PIdatabase_path,
            out_file)
