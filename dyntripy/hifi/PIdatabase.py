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
import pandas as pd
from scipy.signal import spectrogram
from scipy.integrate import simps
import multiprocessing
import warnings

from dyntripy.utils import load_gf, gf


def _psd_integral(pxx_all, f, f_win, gf_parameters):
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


def pi41day(
        sac_file,
        gf_info_file,
        year,
        day,
        sta,
        data_path,
        nperseg,
        noverlap,
        nfft,
        f_win_list,
        out_file):
    try:
        _, sensitivity, normalizing, zeros, poles = load_gf(sac_file, gf_info_file)
        if zeros is not None:
            gf_parameters = [sensitivity, normalizing, zeros, poles]

            st = obspy.read(os.path.join(data_path, str(year), day, '.'.join([day, sta])))
            tr = st[0]
            tr.detrend('linear')
            tr.detrend('constant')
            fs = tr.stats.sampling_rate
            starttime = tr.stats.starttime

            f, t, pxx_all = spectrogram(
                tr.data, fs,
                window='hanning', nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=None,
                return_onesided=True, scaling='density', axis=-1, mode='psd'
            )
            pi_array = np.array([
                [_psd_integral(pxx_all[:, t_i], f, f_win, gf_parameters) for f_win in f_win_list]
                for t_i in range(len(t) - 1)
            ])

            t = np.full(len(t), starttime) + t
            output_array = np.column_stack((t[:-1], pi_array))
            output_df = pd.DataFrame(
                data=output_array,
                columns=['time'] + [str(int(f_win[0])) + '-' + str(int(f_win[1])) for f_win in f_win_list]
            )
            output_df.to_csv(out_file, index=False)

    except Exception as err_msg:
        print(sac_file)
        print(err_msg)
        if os.path.exists(out_file):
            os.remove(out_file)


def pi41day_run(
        sac_file,
        gf_info_file,
        year,
        day,
        sta,
        data_path,
        nperseg,
        noverlap,
        nfft,
        f_win_list,
        out_file
):
    try:
        if os.path.exists(out_file):
            print('%s already exists!' % out_file)
        else:
            if not os.path.exists(os.path.join(data_path, str(year), day, sac_file)):
                warnings.warn('No %s!' % sac_file)
            else:
                pi41day(
                    sac_file,
                    gf_info_file,
                    year,
                    day,
                    sta,
                    data_path,
                    nperseg,
                    noverlap,
                    nfft,
                    f_win_list,
                    out_file
                )
    except Exception as err_msg:
        print(sac_file)
        print(err_msg)
    return None


def run_pi_parallel(
        sta,
        target_dates,
        data_path,
        # time_segment,
        nperseg,
        noverlap,
        nfft,
        f_win_list,
        gf_info_file,
        out_folder,
        cores
):
    pool = multiprocessing.Pool(processes=cores)
    tasks = []

    for target_date in sorted(target_dates):
        year = target_date.year
        day = ''.join([str(target_date.year).zfill(4),
                       str(target_date.month).zfill(2),
                       str(target_date.day).zfill(2)])
        sac_file = '.'.join([day, sta])
        out_file = os.path.join(out_folder, 'PI_' + str(sta) + '_' + str(day) + '.csv')

        # pi41day_run(
        #     sac_file,
        #     gf_info_file,
        #     year,
        #     day,
        #     sta,
        #     data_path,
        #     # time_segment,
        #     nperseg,
        #     noverlap,
        #     nfft,
        #     f_win_list,
        #     out_file
        # )

        tasks.append((
            sac_file,
            gf_info_file,
            year,
            day,
            sta,
            data_path,
            nperseg,
            noverlap,
            nfft,
            f_win_list,
            out_file
        ))

    print('\n')

    # chunksize is how many tasks will be processed by one processor
    rs = pool.starmap_async(pi41day_run, tasks, chunksize=1)
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
