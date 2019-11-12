'''
Written by Yijian Zhou.
The functions:
1. mseed2sac
2. merge
3. change head files
4. slice into days
5. archive (mkdir + mv)

Modified by Naidan YUN on 20180407.

'''
import os, glob
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
import subprocess
from obspy.core import *
from datetime import datetime,timedelta
import shutil
#Tool
def jday2YYYYMMDD(year,jday):
    dt = datetime.strptime(str(year)+str(jday), '%Y%j').date()
    month=str(dt.month)
    day=str(dt.day)
    YYYYMMDD=str(year)+month.zfill(2)+day.zfill(2)

    return YYYYMMDD
def cal_days_for_1year(year):
    dt = datetime.strptime(str(year) + '1231', '%Y%m%d').date()
    jday=int(dt.strftime('%j'))
    return jday
#
def merge(tar_file_list, out_sacfile):
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "r %s\n" % tar_file_list[0]
    print (tar_file_list)
    if len(tar_file_list) > 1:
        for tar_file in tar_file_list[1:]:
            s += "r more %s\n" %tar_file # SAC v101.6 or later
    s += "merge g z o a \n"
    s += "w %s \n" % out_sacfile
    s += "q \n"
    p.communicate(s.encode())
#
def slice_stream_sac(out_sacfile_stage1, tar_ts, tar_te,out_sacfile_stage2):
    str_day=out_sacfile_stage1.split('/')[-2]
    day_date=UTCDateTime(str_day)
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "cuterr fillz \n"
    s += "cut b %s %s \n" %(tar_ts, tar_te)
    s += "r %s \n" % out_sacfile_stage1
    s += "ch b 0 \n"
    s += "ch nzhour %s nzmin %s nzsec %s nzmsec %s \n" %(0,
                                                         0,
                                                         0,
                                                         0)
    s += "ch nzyear %s nzjday %s \n" %(day_date.year, day_date.julday)
    s += "w %s \n" % out_sacfile_stage2
    s += "q \n"
    p.communicate(s.encode())

def slice_getOneday(out_sacfile_stage1,out_sacfile_stage2):
    str_day=out_sacfile_stage1.split('/')[-2]
    day_date=UTCDateTime(str_day)

    st = read(out_sacfile_stage1)
    ts = st[0].stats.starttime
    te = st[0].stats.endtime
    tar_ts=day_date-ts
    tar_te=tar_ts+86400

    slice_stream_sac(out_sacfile_stage1, tar_ts, tar_te,out_sacfile_stage2)
#
def get_tar_day_list(origin_sacfilepath,year,day,sta,chn,firstdate):
    def get_last_sacfile(origin_sacfilepath,year,day,sta,chn):
        last_sacfiles=[]
        for i in range(20):
            today=datetime.strptime(str(year)+str(day), '%Y%j').date()
            tar_date=today-timedelta(days=i+1)
            tar_day=tar_date.strftime('%j')
            tar_year=tar_date.strftime('%Y')
            tar_day_files=glob.glob(os.path.join(origin_sacfilepath,tar_year+'.'+str(tar_day).zfill(3)+'.*.'+sta+'.*.'+chn+'.*'))
            if tar_day_files != []:
                last_sacfiles=last_sacfiles+tar_day_files
                break
        # if last_sacfiles == []:
        #     raise ValueError('Can not find the last target files!')
        return last_sacfiles
    #only need to find the last day
    tar_file_list = []
    today_files = glob.glob(os.path.join(origin_sacfilepath, year + '.' + str(day).zfill(3) + '.*.' + sta + '.*.' + chn + '.*'))

    tar_file_list = tar_file_list + today_files

    if datetime.strptime(str(year)+str(day), '%Y%j').date() != firstdate:
        tar_file_list = tar_file_list + get_last_sacfile(origin_sacfilepath,year,day,sta,chn)

    if tar_file_list == []:
        raise ValueError('There no data for this %s %s'%(year,day))


    return tar_file_list

def merge_and_slice2oneday(origin_sacfilepath,tar_file_list,out_sacfile_stage1,out_sacfile_stage2):
    os.chdir(origin_sacfilepath)
    merge(tar_file_list,out_sacfile_stage1)
    slice_getOneday(out_sacfile_stage1,out_sacfile_stage2)
    os.remove(out_sacfile_stage1)

def run_merge(origin_rootpath, outpath, year_list, sta_list, chn_dic, b_jday_list, e_jday_list,first_flag_list):
    # origin_rootpath='/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Dynamic_triggering_project/Kunming/data/SAC'
    # outpath='/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Dynamic_triggering_project/Kunming/data/raw_new'
    # year_list = \
    #     ['2018']
    # sta_list=['YN.KMI']
    # chn_dic={'YN.KMI':['BHZ','BHE','BHN']}
    # first_jday=269
    # firstdate=datetime.strptime(str(year_list[0]) + str(b_jday_list), '%Y%j').date()
    for y in range(len(year_list)):
        b_jday=b_jday_list[y]
        e_jday=e_jday_list[y]
        if first_flag_list == True:
            firstdate = datetime.strptime(str(year_list[0]) + str(b_jday_list), '%Y%j').date()
        else:
            firstdate=''
        print (firstdate)
        year=year_list[y]
        #year='2009'
        #origin_sacfilepath=os.path.join(origin_rootpath,str(year))
        origin_sacfilepath =origin_rootpath
        out_sacfilepath=os.path.join(outpath,str(year)+'')
        if os.path.exists(out_sacfilepath):
            shutil.rmtree(out_sacfilepath)
        os.makedirs(out_sacfilepath)
        # if y == 0:
        #     day_begin=b_jday_list
        # else:
        #     day_begin = 1
        # days_count = cal_days_for_1year(year)
        # for i in range(day_begin,int(days_count)+1):# the data of one year should be all in the same dir.
        for i in range(int(b_jday), int(e_jday) + 1):
            day=i
            print ('Merge %s %s...'%(str(year),str(day)))
            os.makedirs(os.path.join(out_sacfilepath,jday2YYYYMMDD(str(year),str(day).zfill(3))))

            for sta in sta_list:
                for chn in chn_dic[sta]:
                    #chn='HHZ'
                    tar_file_list = get_tar_day_list(origin_sacfilepath, year,day,sta,chn,firstdate)
                    out_sacfile_stage1=os.path.join(out_sacfilepath,
                                                  jday2YYYYMMDD(str(year),str(day).zfill(3)),
                                                  '.'.join([sta,chn,'pre']))
                    out_sacfile_stage2 = os.path.join(out_sacfilepath,
                                                      jday2YYYYMMDD(str(year), str(day).zfill(3)),
                                                      '.'.join([sta, chn]))
                    merge_and_slice2oneday(origin_sacfilepath,tar_file_list,out_sacfile_stage1,out_sacfile_stage2)

if __name__=="__main__":
    run_merge()
    print ('finish')


