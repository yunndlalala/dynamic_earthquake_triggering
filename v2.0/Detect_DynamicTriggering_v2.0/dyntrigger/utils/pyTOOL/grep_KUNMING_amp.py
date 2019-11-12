'''
Written by Naidan YUN on 20180423.
'''
import os
import re
import pandas as pd
from obspy.core import UTCDateTime
from cal_diatance_YN_sta_QJ import cal_distance

DATA_PATH='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/KUNMING/report/2012-2015'
TELESEISM_CSV_FILE= 'C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/EQ_info.csv'
OUT_FILE= 'C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/match_tele_report_YNsta.csv'

def get_not0(str,end_index):
    i=end_index
    while i >=0 :
        if str[i]==' ':
            begin_index=i+1
            break
        else:
            i=i-1
    out_str=str[begin_index:end_index+1]
    return out_str
def grep_event_phase_amp2csv():
    dis_dic=cal_distance()
    for file in os.listdir(DATA_PATH):
        if file[-4:] == '.txt':
            print ('###############Starting translate %s#############'%file)
            try:
                with open(os.path.join(DATA_PATH,file),'r',encoding='UTF-8') as f:
                    data=f.readlines()
            except:
                with open(os.path.join(DATA_PATH,file),'r',encoding='ansi') as f:
                    data=f.readlines()
            #data_dic={}
            data_array=[]
            event=0
            for line in data:
                if re.split(' +',line)[0]=='DBO':
                    event_time_str=re.split(' +',line)[2]+'T'+re.split(' +',line)[3]+'Z'
                    #event_time_UTC=UTCDateTime(event_time_str)+8*3600
                    #data_dic[event_time_str]=[]
                    event=event_time_str
                    print ('Event is %s.'%event_time_str)
                if re.split(' +',line)[0]=='DPB' and re.split(' +',line)[2] in dis_dic.keys():
                    sta=re.split(' +',line)[2]
                    dis = dis_dic[sta]
                    phase=get_not0(line,27)
                    amp=get_not0(line,86)
                    mag=get_not0(line,104)
                    #data_dic[event].append({phase:amp})
                    data_array.append([event,sta,dis,phase,amp,mag])
            out_dataframe=pd.DataFrame(data=data_array,columns=['event','sta','distance(km)','phase','amp','mag'])
            out_dataframe.to_csv(os.path.join(DATA_PATH,file[:-4]+'.csv'),index=None)
def load_teleseism_data(teleseism_csv_file):
    datafram=pd.read_csv(teleseism_csv_file,encoding='gb2312')
    original_time=[UTCDateTime(time) for time in datafram['time'].values]
    mag=[float(i) for i in datafram['mag'].values]
    result=[[original_time[i],mag[i]] for i in range(len(original_time))]
    return result
def match(teleseism_data):
    out_dataframe_value=[]
    for tele in teleseism_data:
        tar_report_event = [60, None,None, None, 0, None,1000]  # event_diff,report_event,report_sta,report_phase,report_amp,report_mag,distance
        print (tele)
        year=str(tele[0].year)
        report_csv_file=os.path.join(DATA_PATH,year+'.csv')
        report_dataframe = pd.read_csv(report_csv_file)
        for i in range(len(report_dataframe)):
            line=report_dataframe.iloc[i]
            report_event=UTCDateTime(line['event'])
            report_sta = line['sta']
            report_dis=float(line['distance(km)'])
            report_phase = line['phase']
            report_amp = float(line['amp'])
            report_mag=float(line['mag'])
            event_diff=abs(report_event-(tele[0]+3600*8))
            mag_diff=abs(report_mag-tele[1])
            #print (type(report_amp),type(tar_report_event[-2]))
            if event_diff <= tar_report_event[0] and mag_diff <=1:
                if report_sta==tar_report_event[2] and report_amp > tar_report_event[4]:
                    tar_report_event = [event_diff, report_event, report_sta, report_phase, report_amp, report_mag,report_dis]
                if report_sta!=tar_report_event[2] and report_dis < tar_report_event[6]:
                    tar_report_event = [event_diff, report_event, report_sta, report_phase, report_amp, report_mag,report_dis]
        print (tar_report_event)
        out_dataframe_value.append([tele[0],tele[1],tar_report_event[0],tar_report_event[1],
                                    tar_report_event[2],tar_report_event[3],tar_report_event[4],tar_report_event[5],tar_report_event[6]])
    out_dataframe=pd.DataFrame(data=out_dataframe_value,
                               columns=['tele_time','tele_mag','event_diff','report_event','sta','report_phase','report_amp','report_mag','report_dis(km)'])
    out_dataframe.to_csv(OUT_FILE)
if __name__ == '__main__':
    #grep_event_phase_amp2csv()
    teleseism_data=load_teleseism_data(TELESEISM_CSV_FILE)
    match(teleseism_data)
    print ('Finish')


