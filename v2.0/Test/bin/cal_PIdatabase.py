"""
@version:
author:yunnaidan
@time: 2018/10/28
@file: cal_PIdatabase.py
@function:
"""
import os
from dyntrigger.hifi.PIdatabase import run_parallel


def main():
    # Stations
    sta_list = ['NC.GDXB']
    # Channels
    chn_list = ['HHZ']
    # Parameters of the response function.
    sensivity = 1.006000e+09
    normalizing = 9.482280e+18
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
    Gf_parameters = [sensivity, normalizing, zeros, poles]
    # Segment of time for calculating power integral
    time_segment=30 # seconds
    # Frequency windows
    f_win_list = [
        [5, 10],
        [10, 15],
        [15, 20],
        [20, 25],
        [25, 30],
        [30, 35],
        [35, 40]]
    # Path of raw data
    data_path = '../data/raw'
    # Path of result
    out_path = '../result/PI_database'
    # Cores for parallel computing
    cores=1
    for sta in sta_list:
        for chn in chn_list:
            out_folder = os.path.join(out_path, sta + '.' + chn)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            run_parallel(
                sta,
                chn,
                data_path,
                time_segment,
                f_win_list,
                Gf_parameters,
                out_folder,
                cores)


if __name__ == "__main__":
    main()
    print('Finish!')
