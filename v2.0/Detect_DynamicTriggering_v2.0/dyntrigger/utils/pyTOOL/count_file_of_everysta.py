'''
Written by Naidan YUN on 20180415.
In order to count the number of files for every station.
'''
import os
import glob
import datetime
root='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/100HZ'
sta_list=range(1,39)
outfile='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/result/sta_sum.txt'
f=open(outfile,'a')
f.write('sta'+' '+'workday'+' '+'weekend'+'\n')
for sta in sta_list:
    sta_num_weekend=0
    sta_num_workday=0
    sta_name=str(sta).zfill(2)
    files=glob.glob(os.path.join(root,'*','*','*QJ'+sta_name+'*'))
    for file in files:
        day=file.split('\\')[-2]
        date=datetime.datetime.strptime(day,'%Y%m%d')
        if date.isoweekday() in [6,7]:
            sta_num_weekend=sta_num_weekend+1
        else:
            sta_num_workday=sta_num_workday+1
    f.write('QJ'+sta_name+' '+str(sta_num_workday/3)+' '+str(sta_num_weekend/3)+'\n')
f.close()
print ('Finish!!!!!!')