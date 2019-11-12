"""
@version:
author:yunnaidan
@time: 2018/10/29
@file: cal_background_HighPIRation.py
@function:
"""
import os
from dyntrigger.hifi.HighPIRation import run_HighPIRatio


def main():
    # Catalog file of background days
    background_catalog = '../../data/background_catalog.csv'
    # Path of the database saving the results of power integral
    PIdatabase_path = '../../result/PI_database'
    # Path for saving the power integral ratio values of the background days
    out_path='../../result/background_PIR'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # Stations
    sta_list = ['NC.GDXB']
    # Channels
    chn_list=['HHZ']
    for sta in sta_list:
        for chn in chn_list:
            PIdatabase_folder=os.path.join(PIdatabase_path,sta+'.'+chn)
            out_file=os.path.join(out_path,sta+'.'+chn+'_background_PIR.csv')
            run_HighPIRatio(background_catalog, PIdatabase_folder, sta, chn,out_file)

if __name__ == '__main__':
    main()
    print ('Finish!')
