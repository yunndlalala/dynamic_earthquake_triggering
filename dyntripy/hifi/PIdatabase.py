"""
@version:
author:yunnaidan
@time: 2018/10/26
@file: PIdatabase.py
@function:
"""
import os
import time
import obspy
import numpy as np
from scipy.integrate import simps
from obspy.core.utcdatetime import UTCDateTime
import multiprocessing

from dyntripy.utils import psd, load_gf, gf


def psd_integral(pxx_all, f, f_win, gf_parameters):
    f_min = f_win[0]
    f_max = f_win[1]
    f_index = np.where((f >= f_min) & (f <= f_max))
    f_tar = f[f_index]
    pxx_tar = pxx_all[f_index]

    sensitivity, normalizing, zeros, poles = gf_parameters
    gf_list = [gf(sensitivity, normalizing, zeros, poles, f) for f in f_tar]
    counts = len(gf_list)
    pxx_tar_remove_response = [(pxx_tar[i] / (gf_list[i]**2)) * (10**18)
                               for i in range(counts)]  # m/s to nm/s
    pi = simps(pxx_tar_remove_response, f_tar)
    return pi


def abs_time(day_date, time_segment, i):
    year = int(str(day_date)[:4])
    month = int(str(day_date)[4:6])
    day = int(str(day_date)[6:8])
    hour = int((i * time_segment) / 3600)
    minute = int(((i * time_segment) % 3600) / 60)
    second = int(((i * time_segment) % 3600) % 60)
    abs_time_value = UTCDateTime(year, month, day, hour, minute, second)
    return abs_time_value


def pi41day(
        sac_file,
        gf_info_file,
        year,
        day,
        sta,
        data_path,
        time_segment,
        f_win_list,
        out_file):
    try:
        _, sensitivity, normalizing, zeros, poles = load_gf(sac_file, gf_info_file)
        gf_parameters = [sensitivity, normalizing, zeros, poles]

        with open(out_file, 'a') as f:
            f.write('time')
            for f_win in f_win_list:
                f.write(',' + str(int(f_win[0])) + '-' + str(int(f_win[1])))
            f.write('\n')

        st = obspy.read(os.path.join(data_path, str(year),
                                     day, '.'.join([day, sta])))
        tr = st[0]
        tr.detrend('linear')
        tr.detrend('constant')
        fs = tr.stats.sampling_rate
        start_time = tr.stats.starttime

        segment_count = int(86400 / time_segment)
        for i in range(segment_count):
            abs_time_value = abs_time(day, time_segment, i)
            start_point_index = max(0, round((abs_time_value - start_time) * fs))
            end_point_index = max(
                0, round(
                    (abs_time_value + time_segment - start_time) * fs))
            # If the index exceeds the length of the data, no error will be thrown,
            # but empty array.
            data = tr.data[start_point_index:end_point_index]
            if len(data) != 0:
                pxx_all, f = psd(data, fs)
                pi_list = []
                for f_win in f_win_list:
                    pi = psd_integral(pxx_all, f, f_win, gf_parameters)
                    pi_list.append(pi)
            else:
                pi_list = [0.0 for _ in f_win_list]

            with open(out_file, 'a') as f:
                f.write(str(abs_time_value)[:-4] + 'Z')
                for pi in pi_list:
                    f.write(',' + str(pi))
                f.write('\n')
    except Exception as err_msg:
        print(err_msg)


def run_pi_parallel(
        sta,
        target_dates,
        data_path,
        time_segment,
        f_win_list,
        gf_info_file,
        out_folder,
        cores):

    pool = multiprocessing.Pool(processes=cores)
    tasks = []
    for target_date in sorted(target_dates):
        print('Prepare database task ' + str(target_date), end='\r')
        year = target_date.year
        day = ''.join([str(target_date.year).zfill(4),
                       str(target_date.month).zfill(2),
                       str(target_date.day).zfill(2)])
        sac_file = '.'.join([day, sta])
        if not os.path.exists(
            os.path.join(
                data_path,
                str(year),
                day,
                sac_file)):
            continue

        out_file = os.path.join(out_folder, 'PI_' + str(sta) + '_' + str(day) + '.csv')
        if os.path.exists(out_file):
            continue

        # pi41day(year, day, sta, data_path, time_segment, f_win_list, gf_parameters, out_file)
        tasks.append(
            (sac_file,
             gf_info_file,
             year,
             day,
             sta,
             data_path,
             time_segment,
             f_win_list,
             out_file))

    print('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(pi41day, tasks, chunksize=1)
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
