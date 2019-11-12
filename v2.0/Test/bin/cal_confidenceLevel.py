"""
@version:
author:yunnaidan
@time: 2018/10/10
@file: cal_confidenceLevel.py
@function:
"""
import os
from dyntrigger.hifi.ConfidenceLevel import confidence_level_Gauss


def main():
    # Path of the files containing the background PIR after associating
    background_PIR_associated_path = '../result/background_PIR_associated'
    # Path of the files containing the teleseismic PIR
    tele_PIR_path = '../result/tele_PIR'
    # Path for saving the confidence level
    out_path = '../result/confidence_level'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # Stations
    sta_list = ['NC.GDXB']
    # Channels
    chn_list=['HHZ']
    for sta in sta_list:
        for chn in chn_list:
            background_PIR_associated_file=os.path.join(background_PIR_associated_path,sta+'.'+chn+'_background_PIR_associated.csv')
            tele_PIR_file=os.path.join(tele_PIR_path,sta+'.'+chn+'_PIR.csv')
            out_file=os.path.join(out_path,sta+'.'+chn+'_CL.csv')
            confidence_level_Gauss(background_PIR_associated_file, tele_PIR_file, out_file,'PIRatio_log')

    return None


if __name__ == "__main__":
    main()
    print('Finish!')
