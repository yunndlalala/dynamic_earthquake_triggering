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
import logging

logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def cl41event(row, threshold, out_file):
    tele_pir = row['remote_earthquake_ratio']
    if tele_pir != 'None':
        tele_pir = float(tele_pir)
        background_pirs = re.split(' +', row['background_day_ratio'])
        background_pir_list = np.array(background_pirs, dtype='float64')
        if len(background_pir_list) < 30:
            confidence_level_value = 'Insufficient RBs'
            tri_flag = None
        else:
            mean, std = stats.norm.fit(background_pir_list)
            background_norm = stats.norm(loc=mean, scale=std)
            confidence_level_value = background_norm.cdf(tele_pir)
            if confidence_level_value >= threshold:
                tri_flag = int(1)
            else:
                tri_flag = int(0)
    else:
        confidence_level_value = 'No RE'
        tri_flag = None

    with open(out_file, 'a') as f:
        f.write(','.join([str(row['time']), str(confidence_level_value), str(tri_flag)]))
        f.write('\n')

    return None


def run_cl_parallel(
        background_pir_associated_file,
        threshold,
        out_file,
        cores):
    if os.path.exists(out_file):
        os.remove(out_file)
    with open(out_file, 'a') as f:
        f.write(','.join(['time', 'cl', 'triggering']))
        f.write('\n')

    tele_background_pir_df = pd.read_csv(background_pir_associated_file)

    pool = multiprocessing.Pool(processes=cores)
    tasks = []
    for row_index, row in tele_background_pir_df.iterrows():
        print ('Prepare cl task %s' % row['time'], end='\r')
        tasks.append((row, threshold, out_file))
        # cl41event(row, threshold, out_file)
    print ('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(cl41event, tasks, chunksize=1)
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
