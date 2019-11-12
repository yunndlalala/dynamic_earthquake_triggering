'''
Written by Naidan YUN on 20180427.
'''
import os
import obspy

def cal_distance():
    YN_sta_file='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/KUNMING/report/stationlist_YN.txt'
    QJ_lon=102.9990
    QJ_lat=26.7647
    with open(YN_sta_file) as f:
        data=f.readlines()
    lat=[float(line.split(' ')[3]) for line in data]
    lon=[float(line.split(' ')[4]) for line in data]
    sta=[line.split(' ')[1] for line in data]
    dis_dic={}
    for i in range(len(sta)):
        dis=obspy.geodetics.base.gps2dist_azimuth(QJ_lat,QJ_lon,lat[i],lon[i])[0]/1000
        dis_dic[sta[i]]=dis
    return dis_dic
if __name__ == '__main__':
    dis_dic=cal_distance()
    print ('finish')