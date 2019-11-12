"""
@version:
author:yunnaidan
@time: 2018/12/02
@file: visual_varTimewindow.py
@function:
"""
import numpy as np
import pandas as pd
import obspy
import os
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from pylab import *
plt.rc('font',family='Times New Roman')
# plt.rc('xtick', labelsize=10)
# plt.rc('ytick', labelsize=10)
def lastTime(tr):
    i=1
    max_list=[]
    while i*1000 <= len(tr.data):
        data=tr.data[(i-1)*1000:i*1000]
        max_value=np.max(data)
        max_list.append(max_value)
        i=i+1
    max_max_value=np.max(max_list)
    index=max_list.index(max_max_value)
    for i in range(index,len(max_list)):
        if max_list[i]<0.01*max_max_value:
            end_index=index+i
            break
    return (end_index+1)*1000
def lastTime_plot(sacfile,ax):
    st=obspy.read(sacfile)
    tr=st[0]
    fs=tr.stats.sampling_rate
    b=tr.stats.sac['t2']
    b_index=b*fs
    End_index=lastTime(tr)
    # x=np.range(0,len(tr.data),100)
    # y=tr.data[x]
    # ax.plot(x, y, 'k')
    ax.plot(np.arange(0,len(tr.data)-int(b_index),1)/fs,tr.data[int(b_index):],'gray')
    ax.plot([(End_index-b_index)/fs,(End_index-b_index)/fs],[0,np.max(tr.data)],'--',color='gray')

    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
    ax.yaxis.set_major_formatter(xfmt)

    #ax.set_ylabel('Amplitude (nm/s)',fontsize=15)
    #ax.set_xlabel('Time (s)')

    return None
def varTimeWindow_confidenceLevel_plot(confidenceLevel_summaryFile,tele,ax):
    summary_df = pd.read_csv(confidenceLevel_summaryFile)
    x_list = [float(summary_df.columns[i]) for i in range(1,len(summary_df.columns))]
    tar_df = summary_df[summary_df['tele'] == tele]
    y_list = [tar_df[str(int(x))].values[0] for x in x_list]
    ax.plot(x_list, y_list, '-',color='black')
    ax.scatter(x_list, y_list, marker='o', s=100, c='', edgecolors='black')

    min_confidenceLevel=np.nanmin(y_list)
    min_confidenceLevel_index=y_list.index(min_confidenceLevel)
    min_confidenceLevel_time=x_list[min_confidenceLevel_index]
    ax.scatter(min_confidenceLevel_time, min_confidenceLevel, marker='o', s=100, c='black', edgecolors='black')

    text(0.95, 0.9, str(tele),
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax.transAxes,fontsize=15)
    text(0.95, 0.75, 'Minimum confidence level: %.2f'%min_confidenceLevel,
         horizontalalignment='right',
         verticalalignment='center',
         transform=ax.transAxes,fontsize=15)

    #ax.set_ylabel('Significance Level')
    #ax.set_xlabel('Time (s)')
    return None
def oneTimeWindow_confidenceLevel_plot(confidenceLevel_file,catalog_file,tele,ax):

    catalog_df = pd.read_csv(catalog_file)
    B_time = catalog_df[catalog_df['time'] == tele]['B_time'].values[0]
    E_time = catalog_df[catalog_df['time'] == tele]['E_time'].values[0]
    timedelta = UTCDateTime(E_time) - UTCDateTime(B_time)

    tele=str(UTCDateTime(tele))
    confidenceLevel_df=pd.read_csv(confidenceLevel_file)
    tar_df = confidenceLevel_df[confidenceLevel_df['tele'] == tele]
    confidenceLevel=tar_df['confidence_level_value'].values[0]
    ax.scatter(timedelta,confidenceLevel,marker='*',s=1000,c='',edgecolors='black')
    #ax.set_ylabel('Significance Level')
    #ax.set_xlabel('Time (s)')
    return None
if __name__ == '__main__':
    confidenceLevel_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Result/Geysers/vel5_vel2_timewindow/fixed_before_window/25_35frequency_range/norm_distribution/confidenceLevelGDXB.csv'
    catalog_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/vel5_vel2_timewindow/fixed_before_window/25_35frequency_range/EQ_info_GeysersGDXB.csv'
    tele='2010-02-27T06:34:11.53Z'

    date,time=tele.split('T')
    year,month,day=date.split('-')
    hour,min,sec=time[:-1].split(':')
    tele_filename=''.join([year,month,day,hour,min,sec])
    envelope_path='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/removeRES/GDXB/events_short_envelope'
    envelope_file=os.path.join(envelope_path,tele_filename,'*Z')

    confidenceLevel_summaryFile='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Result/Geysers/vel5_varEnd/GDXB/all_confidenceLevelGDXB.csv'

    fig=plt.figure(figsize=[12,8])
    ax=fig.add_subplot(111)

    lastTime_plot(envelope_file, ax)
    ax2 = ax.twinx()
    varTimeWindow_confidenceLevel_plot(confidenceLevel_summaryFile, str(UTCDateTime(tele)), ax2)
    oneTimeWindow_confidenceLevel_plot(confidenceLevel_file, catalog_file, tele, ax2)
    ax.set_xlim([0,18000])
    ax2.set_xlim([0, 18000])
    fig.show()
    print ('Finish')

