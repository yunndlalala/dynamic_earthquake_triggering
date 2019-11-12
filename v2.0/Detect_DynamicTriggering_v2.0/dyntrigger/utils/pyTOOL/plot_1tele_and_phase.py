'''
Written by Naidan YUN on 20180528.
'''
import os
import obspy
import pandas as pd
import numpy as np
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

def main(event_before4hour,phase_file):
    tele_before4hour_datetime=str(UTCDateTime(event_before4hour))[:-5]+'Z'
    phase_dataframe=pd.read_csv(phase_file)
    tar_phase_dataframe=phase_dataframe[phase_dataframe['teleseism']==tele_before4hour_datetime]
    st = obspy.read(os.path.join(TELE_WAVEFORM_FOLDER, event_before4hour, '*BHZ*'))
    st.detrend('linear')
    st.detrend('constant')
    st.filter('bandpass', freqmin=5,freqmax=15)

    fig = plt.figure()
    #ax = fig.add_subplot(1, 1, 1)
    # ax1=fig.add_subplot(111)
    for i in range(len(st)):
        tr=st[i]
        sta = tr.stats.station
        print (sta)
        delta=tr.stats.delta
        phase_time=tar_phase_dataframe[tar_phase_dataframe['sta_name']==sta]['tp'].values
        phaseIndex=phase_time/delta
        wave_data=tr.data
        ax = fig.add_subplot(len(st), 1, i + 1)
        ax.plot([10*i for i in range(int(len(wave_data)/10))], np.array(wave_data[[10*i for i in range(int(len(wave_data)/10))]]), 'black')
        for index in phaseIndex:
            ax.plot([index,index],[np.min(wave_data),np.max(wave_data)],'r--')
        ax.set_yticks([0])
        ax.set_yticklabels([sta])
        ax.set_xticks([500000 * i for i in range((len(wave_data) // 500000) + 1)])
        if i != 0:
            ax.spines['top'].set_color('none')
        if i != len(st) - 1:
            ax.set_xticks([])
            ax.spines['bottom'].set_color('none')
        if i == len(st) - 1:
            ax.set_xticklabels(
                [int((500000 * i) / 100) for i in range((len(wave_data) // 500000) + 1)])
            ax.set_xlabel('Time (s)')
    return fig


if __name__ == "__main__":
    ROOTDIR='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang'
    TELE_WAVEFORM_FOLDER=os.path.join(ROOTDIR,'data/events_v2/events_and_beforeevents_p')
    PHASE_FILE=os.path.join(ROOTDIR,'data/local_phase_time_before_during.csv')
    event_before4hour = '20130715101000.00'
    fig=main(event_before4hour,PHASE_FILE)
    fig.show()
    print ('Finish')


