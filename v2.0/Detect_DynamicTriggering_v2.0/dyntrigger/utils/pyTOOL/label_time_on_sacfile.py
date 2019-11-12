"""
@version:
author:yunnaidan
@time: 2018/10/01
@file: label_time_on_sacfile.py
@function: Label the first motion——a in the head of sac file.
"""
import os
import pandas as pd
import obspy
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
def label_a(data_path,origin_data_file):
    origin_df = pd.read_csv(origin_data_file)
    for index, row in origin_df.iterrows():
        plt.figure()
        time=row['time']
        print(time)
        time_datetime=UTCDateTime(time)

        time_dirname=str(time_datetime.year)+\
                     str(time_datetime.month).zfill(2)+\
                     str(time_datetime.day).zfill(2)+\
                     str(time_datetime.hour).zfill(2)+\
                     str(time_datetime.minute).zfill(2)+\
                     str(time_datetime.second).zfill(2)+\
                     '.'+str(time_datetime.microsecond)[:2]
        st=obspy.read(os.path.join(data_path,time_dirname,'*'))
        delta=st[0].stats.sampling_rate
        depth=row['depth']
        dist_in_km=row['dist(km)']
        distance_in_degree=dist_in_km/111
        arrivals = model.get_travel_times(source_depth_in_km=depth,
                                          distance_in_degree=distance_in_degree)
        travel_time=arrivals[0].time

        for tr in st:
            tr.stats.sac['a']=travel_time
            network=tr.stats.network
            station=tr.stats.station
            channel=tr.stats.channel
            tr.write(os.path.join(data_path,time_dirname,'.'.join([network,station,channel])), format="sac")
        print (travel_time)
        # plt.plot(st[0].data[:100000])
        # plt.plot([travel_time/delta,travel_time/delta],[0,10000])
        # plt.show()
def label_Btime_and_Etime(data_path,origin_data_file):
    origin_df = pd.read_csv(origin_data_file)
    for index, row in origin_df.iterrows():
        time=row['time']
        print(time)
        time_datetime=UTCDateTime(time)

        # time_dirname=str(time_datetime.year)+\
        #              str(time_datetime.month).zfill(2)+\
        #              str(time_datetime.day).zfill(2)+\
        #              str(time_datetime.hour).zfill(2)+\
        #              str(time_datetime.minute).zfill(2)+\
        #              str(time_datetime.second).zfill(2)+\
        #              '.'+str(time_datetime.microsecond)[:3].zfill(3)
        date=time.split('T')[0]
        t=time.split('T')[1]
        time_dirname = ''.join(date.split('-')) + \
                   ''.join(t.split(':'))[:-1]
        st=obspy.read(os.path.join(data_path,time_dirname,'*'))
        A_time = UTCDateTime(row['A_time'])
        B_time=UTCDateTime(row['B_time'])
        E_time=UTCDateTime(row['E_time'])
        for tr in st:
            tr.stats.sac['t1'] = A_time - time_datetime
            tr.stats.sac['t2']=B_time-time_datetime
            tr.stats.sac['t3'] = E_time-time_datetime
            network=tr.stats.network
            station=tr.stats.station
            channel=tr.stats.channel
            tr.write(os.path.join(data_path,time_dirname,'.'.join([network,station,channel])), format="sac")
    return None

if __name__=="__main__":
    origin_data_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/EQ_info.csv'
    data_path='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/event_data/removeRES/events'
    label_Btime_and_Etime(data_path, origin_data_file)
    #label_a(data_path, origin_data_file)
    print ('Finish!')