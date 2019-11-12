"""
@version:
author:yunnaidan
@time: 2018/11/23
@file: associate_background_PIR.py
@function:
"""
import os
from dyntrigger.hifi.AssociateBackgroundPIR import PIR_associate,plot_sort_data,test_windowRatio

def main():
    # Catalog of the telesismic events
    telesesmic_catalog = '../../data/teleseismic_catalog.csv'
    # Path for saving the power integral ratio values of the background days
    background_PIR_path = '../../result/background_PIR'
    # Path for saving the background PIR after associating
    out_path = '../../result/background_PIR_associated'
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    # The number of background days for one teleseismic event
    dayWindow=4 # 120
    # Stations
    sta_list = ['NC.GDXB']
    # Channels
    chn_list=['HHZ']

    for sta in sta_list:
        for chn in chn_list:
            background_PIR_file=os.path.join(background_PIR_path,sta+'.'+chn+'_background_PIR.csv')
            out_file=os.path.join(out_path,sta+'.'+chn+'_background_PIR_associated.csv')
            PIR_associate(telesesmic_catalog, background_PIR_file, 'PIRatio_log', dayWindow, out_file)

if __name__ == '__main__':
    main()
    print ('Finish!')
