import os
import shutil
import glob
import pandas as pd
from datetime import datetime,timedelta
origian_datapath='/home/data/100HZ'
out_datapath='/home/yunnd/Dynamic_triggering/Xiaojiang/data'
catalog='/home/yunnd/Dynamic_triggering/Xiaojiang/data/EQ_info.csv'

def copy_tarfolder(time):
    year = time.year
    month = time.month
    day = time.day
    tar_folder = os.path.join(origian_datapath, str(year), str(year) + '%02d' % month + '%02d' % day)
    
    out_folder = os.path.join(out_datapath, str(year))
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    cmd='scp -r zhouyj@162.105.91.219:%s %s'%(tar_folder,out_folder)
    print (cmd)
    os.system(cmd)

def copy_tarfolder_sta(year,sta_list):
    year_folder = os.path.join(origian_datapath,year)
    for sta in sta_list:
        print ('cp station %s'%sta)
        for day_folder in os.listdir(year_folder):
            tar_files=glob.glob(os.path.join(year_folder,day_folder,'*QJ'+sta+'*'))
            if len(tar_files) != 0:
                out_folder=os.path.join(out_datapath,year,day_folder)
                if not os.path.exists(out_folder):
                    os.mkdir(out_folder)
                for file in tar_files:
                    cmd = 'scp zhouyj@162.105.91.219:%s %s' % (file, out_folder)
                    os.system(cmd)
def run_copy_tarfolder():
    dat=pd.read_csv(catalog)
    for time in dat['time'].values[43:]:
        print (time)
        time=datetime.strptime(time,'%Y-%m-%d %H:%M:%S.%f')
        time2=time-timedelta(days=1)
        time3=time+timedelta(days=1)
        copy_tarfolder(time)
        copy_tarfolder(time2)
        copy_tarfolder(time3)
if __name__ == '__main__':
    copy_tarfolder_sta('2012',['01','','','','','',''])
    print ('Finish!')
