#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2023/05/31
@file: example_data_prepare.py
"""
import os
import time
import numpy as np
import pandas as pd
from glob import glob
import obspy
from dyntripy import utils
import multiprocessing
from datetime import datetime


def rename(file_name):
    net, sta, loc, cha, year, jday = file_name.split('.')
    date = datetime.strptime(year + '.' + jday, '%Y.%j')
    new_date = ''.join([
        str(date.year).zfill(4),
        str(date.month).zfill(2),
        str(date.day).zfill(2),
    ])
    new_name = '.'.join([new_date, net, sta, loc, cha])
    return new_date, new_name


def merge_1day(
        file,
        merged_data_path,
        merged_mseed_file
):
    """
    Merge the data of one day and save it to the merged_data_path, if the
    data is several segments of one day.
    """
    try:
        new_date, new_name = rename(file.split('/')[-1])
        merged_data_folder = os.path.join(merged_data_path, new_date)
        if not os.path.exists(merged_data_folder):
            os.makedirs(merged_data_folder)
        output_file = os.path.join(merged_data_folder, new_name)

        if not os.path.exists(output_file):

            st = obspy.read(file)

            if len(st) > 1:

                with open(merged_mseed_file, 'a+') as f:
                    f.write(file + '\n')
                st.merge(method=1, fill_value=0)

            st.write(output_file, format='MSEED')

            os.remove(file)

    except Exception as err_msg:
        print(err_msg)
    return None


def merge_1day_parallel():
    """
    Same as merge_1day, but use multiprocessing to speed up the process.
    """
    raw_data_path = '/data/DATA'
    merged_data_path = '/data/DATA'
    merged_mseed_file = '/data/fault_mseed_BHEN.txt'
    if not os.path.exists(merged_data_path):
        os.makedirs(merged_data_path)

    cores = 48
    pool = multiprocessing.Pool(processes=cores)
    tasks = []

    files = sorted(glob(os.path.join(raw_data_path, '*/*')))
    for file in files:
        if file[0] != '.' and len(file.split('.')) == 6:
            # merge_1day(file, merged_data_path, merged_mseed_file)
            tasks.append((file, merged_data_path, merged_mseed_file))

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(merge_1day, tasks, chunksize=1)
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

    return None


def catalog_processing():
    """
    Process the raw catalog file to get new catalog file as the input.
    """
    output = '../data/remote_earthquake_catalog.csv'
    catalog_file = '../data/raw_remote_earthquake_catalog.csv'
    sta_file = '../data/station.csv'
    utils.gen_time_windows(
        catalog_file=catalog_file,  # the raw catalog file
        sta_file=sta_file,  # the station file
        tb=18000,  # the time window before the remote earthquake, unit: second
        te_b_vel=5.0,  # the velocity of a surface wave whose arrival is the start of $T_e$
        te_e_vel=2.0,  # the velocity of a surface wave whose arrival is the end of $T_e$
        f_min=5,  # the minimum frequency of the high-frequency band
        f_max=20,  # the maximum frequency of the high-frequency band
        f_bin=5,  # the frequency bin
        out_file=output,  # the output file
    )
    return None


def gen_pz_summary():
    """
    Generate the summary of the PZ files and the waveform files as the input.
    """
    def catch_tr_name(pz_file):
        _, _, net, sta, cha, loc, _, _ = pz_file.split('_')
        return '.'.join([net, sta, loc, cha])

    waveform_path = '/data/DATA'
    pz_root_path = '../data/pz/rdseed_pz'
    output_file = '../data/sacPZ.csv'

    pz_files = glob(os.path.join(pz_root_path, '*'))
    pz_tr_list = np.array([catch_tr_name(file.split('/')[-1]) for file in pz_files])

    sac_info = []
    for root, dirs, files in sorted(os.walk(waveform_path)):
        for file in sorted(files):
            if file[0] != '.':
                print(file)
                date, net, sta, loc, cha = file.split('.')
                tr_name = '.'.join([net, sta, loc, cha])
                pz_index = np.where(pz_tr_list == tr_name)[0]
                if len(pz_index) == 0:
                    pz_file = None
                else:
                    pz_file = pz_files[pz_index[0]]
                sac_info.append([
                    file, pz_file
                ])

    out_df = pd.DataFrame(data=sac_info, columns=['sacfile', 'PZ_file'])

    out_df.to_csv(output_file, index=False)

    return None


if __name__ == '__main__':

    # merge_1day_parallel()

    # catalog_processing()

    pass