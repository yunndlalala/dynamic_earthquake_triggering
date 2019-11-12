"""
@version:
author:yunnaidan
@time: 2018/10/29
@file: HighPIRation.py
@function:
"""
import os
import numpy as np
import pandas as pd
import time
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
from obspy.core.utcdatetime import UTCDateTime
import logging
logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


# def load_PIdata(sta, t1, t4, data_path):
#     B_day = ''.join([str(t1.year),
#                      str(t1.month).zfill(2),
#                      str(t1.day).zfill(2)])
#     E_day = ''.join([str(t4.year),
#                      str(t4.month).zfill(2),
#                      str(t4.day).zfill(2)])
#     file_exist = False
#     if B_day == E_day:
#         tar_file = os.path.join(
#             data_path,
#             str(sta) +
#             '_' +
#             str(B_day) +
#             '.csv')
#         if os.path.exists(tar_file):
#             file_exist = True
#             PIdata_frame = pd.read_csv(
#                 os.path.join(
#                     data_path,
#                     str(sta) +
#                     '_' +
#                     str(B_day) +
#                     '.csv'))
#         else:
#             PIdata_frame = None
#     else:
#         tar_file1 = os.path.join(
#             data_path, str(sta) + '_' + str(B_day) + '.csv')
#         tar_file2 = os.path.join(
#             data_path, str(sta) + '_' + str(E_day) + '.csv')
#         if os.path.exists(tar_file1) and os.path.exists(tar_file2):
#             file_exist = True
#             PIdata_frame1 = pd.read_csv(
#                 os.path.join(
#                     data_path,
#                     str(sta) +
#                     '_' +
#                     str(B_day) +
#                     '.csv'))
#             PIdata_frame2 = pd.read_csv(
#                 os.path.join(
#                     data_path,
#                     str(sta) +
#                     '_' +
#                     str(E_day) +
#                     '.csv'))
#             PIdata_frame = PIdata_frame1.append(
#                 PIdata_frame2, ignore_index=True)
#         else:
#             PIdata_frame = None
#     return file_exist, PIdata_frame

def load_PIdata(sta, chn, Tb_B, Te_E, PIdatabase_path):
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
            str(sta) +'.'+chn+
            '_' +
            str(B_day) +
            '.csv')
        if os.path.exists(tar_file):
            file_exist = True
            PIdata_frame = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) +'.'+chn+
                    '_' +
                    str(B_day) +
                    '.csv'))
        else:
            PIdata_frame = None
    else:
        tar_file1 = os.path.join(
            PIdatabase_path, str(sta) + '.'+chn+'_' + str(B_day) + '.csv')
        tar_file2 = os.path.join(
            PIdatabase_path, str(sta) + '.'+chn+'_' + str(E_day) + '.csv')
        if os.path.exists(tar_file1) and os.path.exists(tar_file2):
            file_exist = True
            PIdata_frame1 = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) +'.'+chn+
                    '_' +
                    str(B_day) +
                    '.csv'))
            PIdata_frame2 = pd.read_csv(
                os.path.join(
                    PIdatabase_path,
                    str(sta) +'.'+chn+
                    '_' +
                    str(E_day) +
                    '.csv'))
            PIdata_frame = PIdata_frame1.append(
                PIdata_frame2, ignore_index=True)
        else:
            PIdata_frame = None
    return file_exist, PIdata_frame


def get_t1234_timedelta(A_time, B_time, E_time):
    t3 = B_time
    t3_int = UTCDateTime(t3.year, t3.month, t3.day, t3.hour, t3.minute)
    t4 = E_time
    win = t4 - t3_int

    t2 = A_time
    t2_int = UTCDateTime(t2.year, t2.month, t2.day, t2.hour, t2.minute)
    t1 = t2_int + 60 - 3600 * 5
    return t1, t2_int, t3_int, t4


def f_column_name(f_min, f_max):
    if f_min == 5 and f_max == 15:
        column_name = ['5-10', '10-15']
    if f_min == 10 and f_max == 20:
        column_name = ['10-15', '15-20']
    if f_min == 15 and f_max == 25:
        column_name = ['15-20', '20-25']
    if f_min == 20 and f_max == 30:
        column_name = ['20-25', '25-30']
    if f_min == 25 and f_max == 35:
        column_name = ['25-30', '30-35']
    if f_min == 30 and f_max == 40:
        column_name = ['30-35', '35-40']
    if f_min == 5 and f_max == 35:
        column_name = ['5-10','10-15','15-20','20-25','25-30','30-35']
    if f_min == 10 and f_max == 35:
        column_name = ['10-15','15-20','20-25','25-30','30-35']
    if f_min == 15 and f_max == 35:
        column_name = ['15-20','20-25','25-30','30-35']
    if f_min == 20 and f_max == 35:
        column_name = ['20-25','25-30','30-35']
    return column_name



def sum_PI(column_name, tar_df):
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
            # This modify has been applied and tested, the changes of results are really small.
            PI_mean = np.mean(PI_list)
            PI_mean_list.append(PI_mean)
    if interrupt:
        PI = 0
    else:
        PI = np.sum(PI_mean_list)


    list4plot=np.zeros(len(tar_df))
    for col in column_name:
        PI_list = [float(i) for i in tar_df[col].values]
        list4plot += PI_list

    return interrupt, PI,list4plot


def PIRatio(PI_b, PI_e):
    # The modify in function sum_PI makes this time normalization not needed.
    # Noted that the time normalization commented out here is not very right.
    # PI_e_2=PI_e/(t4-t3_int)
    # PI_b_2 = PI_b / (t2_int+60 - t1)
    PI_e_2 = PI_e
    PI_b_2 = PI_b
    return PI_e_2 / PI_b_2


def HighPIRatio_for_1event(
        sta,
        chn,
        event,
        A_time,
        B_time,
        E_time,
        f_min,
        f_max,
        data_path,
        outfile):
    t1, t2_int, t3_int, t4 = get_t1234_timedelta(A_time, B_time, E_time)
    file_exist, PIdata_frame = load_PIdata(sta, chn,t1, t4, data_path)
    if file_exist:
        PIdata_frame['Time'] = [UTCDateTime(Time)
                                for Time in PIdata_frame['Time']]
        tar_df_b = PIdata_frame[(PIdata_frame['Time'] > t1 - 60)
                                & (PIdata_frame['Time'] <= t2_int)]
        tar_df_e = PIdata_frame[(PIdata_frame['Time'] >= t3_int) & (
            PIdata_frame['Time'] < t4)]

        column_name = f_column_name(f_min, f_max)
        interrupt_b, PI_b,list4plot_b = sum_PI(column_name, tar_df_b)
        interrupt_e, PI_e,list4plot_e = sum_PI(column_name, tar_df_e)

        str_list4plot_b=';'.join([str(i) for i in list4plot_b])
        str_list4plot_e = ';'.join([str(i) for i in list4plot_e])
        if interrupt_b or interrupt_e:
            interrupt_final = True
            PI_ratio = 0
        else:
            interrupt_final = False
            # PI_ratio = PIRatio(PI_b, PI_e)
            PI_ratio = PIRatio(PI_b, PI_e)
            PI_ratio_log=np.log10(PI_ratio)

        if interrupt_final:
            logging.warning(
                'There is interrupt in this event %s!' %
                str(event))
        else:
            with open(outfile, 'a') as f:
                # f.write(
                #     ','.join([str(event), str(PI_b), str(PI_e), str(PI_ratio),str_list4plot_b,str_list4plot_e]))
                f.write(
                    ','.join([str(event), str(PI_b), str(PI_e), str(PI_ratio),str(PI_ratio_log)]))
                f.write('\n')

        #plt.plot(list4plot_b,'-b')
        #plt.plot(range(len(list4plot_b),len(list4plot_b)+len(list4plot_e)),list4plot_e,'-r')


    # If there is no data for the day from Tb to Te, we will through this
    # distant event.
    else:
        logging.warning('There is no data for this event %s!' % str(event))




def run_HighPIRatio(catalog, data_path, sta, chn,outfile):
    if os.path.exists(outfile):
        os.remove(outfile)
    with open(outfile, 'a') as f:
        f.write(','.join(['event', 'PI_b', 'PI_e', 'PIRatio','PIRatio_log']))
        f.write('\n')

    event_dataframe = pd.read_csv(catalog)
    length = len(event_dataframe)
    for i in range(length):
        #i=45
        event = UTCDateTime(event_dataframe.iloc[i]['time'])
        print('Calculate %s...' % event)
        A_time = UTCDateTime(event_dataframe.iloc[i]['A_time'])
        B_time = UTCDateTime(event_dataframe.iloc[i]['Te_B_time'])
        E_time = UTCDateTime(event_dataframe.iloc[i]['Te_E_time'])
        f_min = event_dataframe.iloc[i]['f_min']
        f_max = event_dataframe.iloc[i]['f_max']
        HighPIRatio_for_1event(
            sta,
            chn,
            event,
            A_time,
            B_time,
            E_time,
            f_min,
            f_max,
            data_path,
            outfile)



if __name__ == '__main__':
    catalog = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Data/Geysers/EQ_info_GeysersGDXB.csv'
    data_path = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Result/Geysers/Energy_Ratio/PIdatabase/all_day'
    sta = 'NC.GDXB'
    outfile = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Result/Geysers/Energy_Ratio/tele_HighPID_GDXB.csv'
    run_HighPIRatio(catalog, data_path, sta, outfile)
