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
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.integrate import simps
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
from obspy.core.utcdatetime import UTCDateTime
import multiprocessing
import logging

logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')

# Plot response function


def plot_Gf():
    sensivity = 1
    zeros = [0, 0, 0, 0]
    poles = [-5.0265 + 3.76987j,
             -5.0265 - 3.76987j,
             -0.5969,
             -0.5969,
             -276.46,
             -276.46,
             -48.0915 + 116.097j,
             -48.0915 - 116.097j,
             -116.101 + 48.0832j,
             -116.101 - 48.0832j]
    normalizing = 1.936754 * (10 ** 13)
    f_list = np.arange(0.001, 1000, 0.001)
    Gf_list = []
    for f in f_list:
        Gf_list.append(Gf(sensivity, normalizing, zeros, poles, f))
    print(Gf(sensivity, normalizing, zeros, poles, 1))
    plt.plot(f_list, Gf_list)
    plt.plot([0.001, 100], [1, 1])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


def psd(data, fs):
    f,pxx_all=welch(data,
                     fs,
                     window='hanning',
                     nperseg=512,
                     noverlap=256,
                     nfft=512,
                     detrend=None,
                     return_onesided=True,
                     scaling='density',
                     axis=-1)
    return pxx_all, f


def Gf(sensivity, normalizing, zeros, poles, f):
    s = complex(0, 2 * np.pi * f)
    Gf1 = 1
    for zero in zeros:
        Gf1 = Gf1 * (s - zero)
    Gf2 = 1
    for pole in poles:
        Gf2 = Gf2 * (s - pole)
    Gf = sensivity * normalizing * Gf1 / Gf2

    return abs(Gf)


def psd_interal(pxx_all, f, f_win, Gf_parameters):
    f_min = f_win[0]
    f_max = f_win[1]
    f_index = np.where((f >= f_min) & (f <= f_max))
    f_tar = f[f_index]
    pxx_tar = pxx_all[f_index]

    sensivity, normalizing, zeros, poles = Gf_parameters
    Gf_list = [Gf(sensivity, normalizing, zeros, poles, f) for f in f_tar]
    counts = len(Gf_list)
    pxx_tar_removeGf = [(pxx_tar[i] / (Gf_list[i]**2)) * (10**18)
                        for i in range(counts)]  # m/s to nm/s
    PI = simps(pxx_tar_removeGf, f_tar)
    return PI


def absTime(day_date, time_segment,i):
    year = int(str(day_date)[:4])
    month = int(str(day_date)[4:6])
    day = int(str(day_date)[6:8])
    hour = int((i*time_segment)/ 3600)
    minute = int(((i*time_segment) % 3600)/60)
    second=int(((i*time_segment) % 3600)%60)
    absTime = UTCDateTime(year, month, day, hour, minute,second)
    return absTime


def run_cal(
        year,
        day,
        sta,
        chn,
        data_path,
        time_segment,
        f_win_list,
        Gf_parameters,
        out_folder):

    outfile = os.path.join(out_folder, str(sta) +'.'+ str(chn)+'_' + str(day) + '.csv')

    if os.path.exists(outfile):
        os.remove(outfile)
    with open(outfile, 'a') as f:
        f.write('Time')
        for f_win in f_win_list:
            f.write(',' + str(f_win[0]) + '-' + str(f_win[1]))
        f.write('\n')

    st = obspy.read(os.path.join(data_path, year, day, '.'.join([sta, chn])))
    tr = st[0]
    tr.detrend('linear')
    tr.detrend('constant')
    fs = tr.stats.sampling_rate

    segment_count=int(86400/time_segment)
    for i in range(segment_count):
        abs_time = absTime(day, time_segment,i)
        data = tr.data[int(i * time_segment * fs):int((i + 1) * time_segment * fs + 1)]
        pxx_all, f = psd(data, fs)
        PI_list = []
        for f_win in f_win_list:
            PI = psd_interal(pxx_all, f, f_win, Gf_parameters)
            PI_list.append(PI)

        with open(outfile, 'a') as f:
            #f.write(str(day) + ',' + str(i))
            f.write(str(abs_time))
            for PI in PI_list:
                f.write(',' + str(PI))
            f.write('\n')


def run_parallel(sta, chn, data_path, time_segment,f_win_list, Gf_parameters, out_folder, cores):

    for year in os.listdir(data_path):
        print('Calculate %s.%s in %s' % (str(sta), str(chn),str(year)))
        pool = multiprocessing.Pool(processes=cores)
        tasks = []
        for day in os.listdir(os.path.join(data_path, year)):
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
        rs = pool.starmap_async(run_cal, tasks, chunksize=1)
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


if __name__ == "__main__":
    plot_Gf()
    sta = 'NC.GDXB'
    chn = 'HHZ'
    data_path = '/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data_local/Geysers/RAW_data/allData_test'
    f_win_list = [[5, 10], [10, 15], [15, 20], [20, 25]]
    sensivity = 1.006000e+09
    zeros = [0,
             0,
             complex(-4.618140e+02, -4.290780e+02),
             complex(-4.618140e+02, 4.290780e+02),
             complex(-1.946970e+02, 0),
             complex(-1.514870e+01, 0)]
    poles = [complex(-3.681800e-02, -3.692300e-02),
             complex(-3.681800e-02, 3.692300e-02),
             complex(-1.023970e+04, -2.725020e+03),
             complex(-1.023970e+04, 2.725020e+03),
             complex(-9.512740e+03, -1.146990e+04),
             complex(-9.512740e+03, 1.146990e+04),
             complex(-4.545250e+02, 0.000000e+00),
             complex(-3.981840e+02, 0.000000e+00),
             complex(-9.294710e+01, -3.889920e+02),
             complex(-9.294710e+01, 3.889920e+02),
             complex(-1.545030e+01, 0.000000e+00)]
    #zeros = [zero / (2 * np.pi) for zero in poles]
    #poles=[pole/(2*np.pi) for pole in poles]
    # n=len(poles)-len(zeros)
    normalizing = 9.482280e+18
    # normalizing=normalizing/((2*np.pi)**n)

    outpath = '/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Result/Geysers/Energy_Ratio/PIdatabase'
    Gf_parameters = [sensivity, normalizing, zeros, poles]
    run_parallel(sta, chn, data_path, f_win_list, Gf_parameters, outpath)
    print('finish')
