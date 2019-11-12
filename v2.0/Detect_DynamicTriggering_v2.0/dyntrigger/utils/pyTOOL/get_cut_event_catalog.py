'''
Created by Naidan YUN at 2018/03/13
The script is in order to grep the basic information for cut_event.py from EQ_info.csv
'''
import os
import pandas as pd
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

def get_catalog(origin_data_file,out_file,time_win_before=-100,time_win_after=100,cal_Parrivetime=True):
    origin_df = pd.read_csv(origin_data_file,encoding="gb2312")
    out_df_data=[]
    for index, row in origin_df.iterrows():
        time=row['time']
        print(time)
        time_datetime=UTCDateTime(time)
        depth=row['depth']
        dist_in_km=row['dist(km)']
        distance_in_degree=dist_in_km/111

        if cal_Parrivetime:
            arrivals = model.get_travel_times(source_depth_in_km=depth,
                                              distance_in_degree=distance_in_degree)
            print (arrivals[0])
            time_win_before_final=time_win_before+arrivals[0].time
            time_win_after_final = time_win_after + arrivals[0].time
        out_df_data.append([time,row['latitude'],row['longitude'],depth,row['mag'],time_win_before_final,time_win_after_final])
    out_df=pd.DataFrame(data=out_df_data,columns=['time','latitude','longitude','depth','mag','time_before','time_after'])
    out_df.to_csv(out_file,index=False)


if __name__=="__main__":
    origin_data_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/EQ_info.csv'
    out_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/catalog4cut_event.csv'
    get_catalog(origin_data_file, out_file, time_win_before=-14400,time_win_after=14400,cal_Parrivetime=True)
    print ('Finish!')