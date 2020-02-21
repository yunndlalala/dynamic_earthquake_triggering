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
import logging

from dyntrigger.utils.basic_utils import psd, load_gf, gf


logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


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
    abs_time = UTCDateTime(year, month, day, hour, minute, second)
    return abs_time


def pi41day(
        year,
        day,
        sta,
        chn,
        data_path,
        time_segment,
        f_win_list,
        gf_parameters,
        out_folder):

    outfile = os.path.join(
        out_folder,
        str(sta) +
        '.' +
        str(chn) +
        '_' +
        str(day) +
        '.csv')
    with open(outfile, 'a') as f:
        f.write('time')
        for f_win in f_win_list:
            f.write(',' + str(int(f_win[0])) + '-' + str(int(f_win[1])))
        f.write('\n')

    st = obspy.read(os.path.join(data_path, year,
                                 day, '.'.join([day, sta, chn])))
    tr = st[0]
    tr.detrend('linear')
    tr.detrend('constant')
    fs = tr.stats.sampling_rate

    segment_count = int(86400 / time_segment)
    for i in range(segment_count):
        abs_time_value = abs_time(day, time_segment, i)
        # If the index exceeds the length of the data, no error will be thrown,
        # but empty array.
        data = tr.data[int(i * time_segment * fs):int((i + 1) * time_segment * fs + 1)]
        pxx_all, f = psd(data, fs)
        pi_list = []
        for f_win in f_win_list:
            pi = psd_integral(pxx_all, f, f_win, gf_parameters)
            pi_list.append(pi)

        with open(outfile, 'a') as f:
            f.write(str(abs_time_value))
            for pi in pi_list:
                f.write(',' + str(pi))
            f.write('\n')


def run_pi_parallel(
        sta,
        chn,
        data_path,
        time_segment,
        f_win_list,
        gf_info_file,
        out_folder,
        cores):

    for year in sorted(os.listdir(data_path)):
        print('Calculate %s.%s in %s' % (str(sta), str(chn), str(year)))
        pool = multiprocessing.Pool(processes=cores)
        tasks = []
        for day in sorted(os.listdir(os.path.join(data_path, year))):
            print(day, end='\r')
            sac_file = '.'.join([day, sta, chn])
            if not os.path.exists(os.path.join(data_path, year, day, sac_file)):
                continue

            _, sensitivity, normalizing, zeros, poles = load_gf(sac_file, gf_info_file)
            gf_parameters = [sensitivity, normalizing, zeros, poles]
            # pi41day(year, day, sta, chn, data_path, time_segment, f_win_list, gf_parameters, out_folder)
            tasks.append(
                (year,
                 day,
                 sta,
                 chn,
                 data_path,
                 time_segment,
                 f_win_list,
                 gf_parameters,
                 out_folder))

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
