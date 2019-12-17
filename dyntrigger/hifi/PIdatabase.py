"""
@version:
author:yunnaidan
@time: 2018/10/26
@file: PIdatabase.py
@function:
"""
import os
import re
import time
import obspy
import math
import pandas as pd
import numpy as np
from scipy.signal import welch
from scipy.integrate import simps
from obspy.core.utcdatetime import UTCDateTime
import multiprocessing
import logging


logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def psd(data, fs):
    nfft = 512
    seg = math.ceil(len(data) / (nfft / 2))
    nperseg = int(len(data) / seg) * 2
    f, pxx_all = welch(data,
                       fs,
                       window='hanning',
                       nperseg=nperseg,
                       noverlap=nperseg / 2,
                       nfft=nfft,
                       detrend=None,
                       return_onesided=True,
                       scaling='density',
                       axis=-1)

    return pxx_all, f


def gf(sensivity, normalizing, zeros, poles, f):
    s = complex(0, 2 * np.pi * f)
    Gf1 = 1
    for zero in zeros:
        Gf1 = Gf1 * (s - zero)
    Gf2 = 1
    for pole in poles:
        Gf2 = Gf2 * (s - pole)
    Gf = sensivity * normalizing * Gf1 / Gf2

    return abs(Gf)


def psd_integral(pxx_all, f, f_win, Gf_parameters):
    f_min = f_win[0]
    f_max = f_win[1]
    f_index = np.where((f >= f_min) & (f <= f_max))
    f_tar = f[f_index]
    pxx_tar = pxx_all[f_index]

    sensivity, normalizing, zeros, poles = Gf_parameters
    Gf_list = [gf(sensivity, normalizing, zeros, poles, f) for f in f_tar]
    counts = len(Gf_list)
    pxx_tar_removeGf = [(pxx_tar[i] / (Gf_list[i]**2)) * (10**18)
                        for i in range(counts)]  # m/s to nm/s
    PI = simps(pxx_tar_removeGf, f_tar)
    return PI


def abs_time(day_date, time_segment, i):
    year = int(str(day_date)[:4])
    month = int(str(day_date)[4:6])
    day = int(str(day_date)[6:8])
    hour = int((i * time_segment) / 3600)
    minute = int(((i * time_segment) % 3600) / 60)
    second = int(((i * time_segment) % 3600) % 60)
    absTime = UTCDateTime(year, month, day, hour, minute, second)
    return absTime


def load_gf(sacfile, Gf_info_file):
    df = pd.read_csv(Gf_info_file)
    PZ_file = df.loc[(df['sacfile'].values == sacfile), 'PZ_file'].values[0]

    with open(PZ_file, 'r') as f:
        text = f.readlines()

    sensivity_line = text[21]
    sensivity = float(sensivity_line.split(':')[1].split('(')[0])
    normalizing_line = text[22]
    normalizing = float(normalizing_line.split(':')[1][:-1])

    zeroNo = int(re.split(' +', text[24])[1])
    zeros = []
    for i in range(zeroNo):
        zero_info = text[25 + i]
        _, real, im = re.split(' +', zero_info[:-1])
        zeros.append(complex(float(real), float(im)))

    pole_line_index = 24 + zeroNo + 1
    poleNo = int(re.split(' +', text[pole_line_index])[1])
    poles = []
    for i in range(poleNo):
        pole_info = text[pole_line_index + 1 + i]
        _, real, im = re.split(' +', pole_info[:-1])
        poles.append(complex(float(real), float(im)))

    Gf_parameters = [sensivity, normalizing, zeros, poles]

    return Gf_parameters


def pi41day(
        year,
        day,
        sta,
        chn,
        data_path,
        time_segment,
        f_win_list,
        Gf_parameters,
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
            f.write(',' + str(f_win[0]) + '-' + str(f_win[1]))
        f.write('\n')

    st = obspy.read(os.path.join(data_path, year,
                                 day, '.'.join([day, sta, chn])))
    tr = st[0]
    tr.detrend('linear')
    tr.detrend('constant')
    fs = tr.stats.sampling_rate

    segment_count = int(86400 / time_segment)
    for i in range(segment_count):
        abs_time = abs_time(day, time_segment, i)
        # If the index exceeds the length of the data, no error will be thrown,
        # but empty array.
        data = tr.data[int(i * time_segment * fs):int((i + 1) * time_segment * fs + 1)]
        pxx_all, f = psd(data, fs)
        PI_list = []
        for f_win in f_win_list:
            PI = psd_integral(pxx_all, f, f_win, Gf_parameters)
            PI_list.append(PI)

        with open(outfile, 'a') as f:
            f.write(str(abs_time))
            for PI in PI_list:
                f.write(',' + str(PI))
            f.write('\n')


def run_pi_parallel(
        sta,
        chn,
        data_path,
        time_segment,
        f_win_list,
        Gf_info_file,
        out_folder,
        cores):

    for year in sorted(os.listdir(data_path)):
        print('Calculate %s.%s in %s' % (str(sta), str(chn), str(year)))
        pool = multiprocessing.Pool(processes=cores)
        tasks = []
        for day in sorted(os.listdir(os.path.join(data_path, year))):
            print(day, end='\r')
            sacfile = '.'.join([day, sta, chn])
            if not os.path.exists(os.path.join(data_path, year, day, sacfile)):
                continue

            Gf_parameters = load_gf(sacfile, Gf_info_file)
            tasks.append(
                (year,
                 day,
                 sta,
                 chn,
                 data_path,
                 time_segment,
                 f_win_list,
                 Gf_parameters,
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
