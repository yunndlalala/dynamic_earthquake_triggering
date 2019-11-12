'''
Written by Naidan YUN on 20180407.
Modified by Naidan YUN on 20181006.

'''
import obspy
import numpy as np
import itertools
import time
import os
import shutil
import glob
import logging
import multiprocessing

# Interrupt checking
def checkInterrupt(trace,threshold_in_point=6000):
    def interrupt_len(data,threshold_in_point,amp_thre=0.00001):
        sum_list=[]
        sum_list_append = sum_list.append
        #data_length=len(data)
        flag=False
        for item in data:
            if abs(item) <= amp_thre:
                sum_list_append(item)
                if len(sum_list) >= threshold_in_point:
                    flag = True
                    break
            else:
                sum_list = []
                sum_list_append = sum_list.append

        return flag

    data=trace.data
    result=interrupt_len(data,threshold_in_point,amp_thre=0.00001)
    return result
def checkInterrupt_deleteOneday(year_folder,out_year_folder,threshold_in_sec=60):
    def sum_sta(day_folder):
        sta_list=[]
        for sac_file in os.listdir(day_folder):
            sta=sac_file[:-4]
            if sta not in sta_list:
                sta_list.append(sta)
        return sta_list

    if not os.path.exists(out_year_folder):
        os.makedirs(out_year_folder)

    for day_folder in sorted(os.listdir(year_folder)):
        sta_list=sum_sta(os.path.join(year_folder,day_folder))

        delete_sta_list=[]
        for sta in sta_list:
            interrupt = False
            st=obspy.read(os.path.join(year_folder,day_folder,'*'+sta+'*'))
            if checkInterrupt(st[0],threshold_in_sec=threshold_in_sec) \
                    or checkInterrupt(st[1],threshold_in_sec=threshold_in_sec) \
                    or checkInterrupt(st[2],threshold_in_sec=threshold_in_sec):
                interrupt = True
            if interrupt:
                delete_sta_list.append(sta)
                if not os.path.exists(os.path.join(out_year_folder,day_folder)):
                    os.makedirs(os.path.join(out_year_folder,day_folder))
                for file in glob.glob(os.path.join(year_folder,day_folder,'*'+sta+'*')):
                    shutil.move(file,os.path.join(out_year_folder,day_folder))
        logging.info('delete the data of these stations on %s: %s'%(day_folder,str(delete_sta_list)))
def checkLack_Interrupt_recordOneday(year_folder,result_file,threshold_in_point=6000):
    def sum_sta(day_folder):
        sta_list=[sac_file[:-4] for sac_file in os.listdir(day_folder)]
        return list(set(sta_list))

    for day_folder in os.listdir(year_folder):
        print (day_folder)
        sta_list=sum_sta(os.path.join(year_folder,day_folder))

        delete_sta_list=[]
        delete_sta_list_append=delete_sta_list.append
        for sta in sta_list:
            s = time.time()
            files = glob.glob(os.path.join(year_folder, day_folder,''.join(['*',sta,'*'])))
            if len(files) < 3:
                delete_sta_list_append(sta)
                continue
            e = time.time()
            print ('pre %f'%(e-s))
            s = time.time()
            st=obspy.read(os.path.join(year_folder,day_folder,''.join(['*',sta,'*'])))
            e = time.time()
            print('read %f'%(e - s))
            s = time.time()
            for tr in st:
                if checkInterrupt(tr,threshold_in_point=threshold_in_point):
                    delete_sta_list_append(sta)
                    break
            e = time.time()
            print('cal %f'%(e - s))
        s = time.time()
        if delete_sta_list != []:
            with open(result_file, 'a') as f:
                f.write(' '.join([str(year_folder),str(day_folder),str(delete_sta_list),'\n']))
        e = time.time()
        print('write %f'%(e - s))
def checkInterrupt_day(day_folder,year_folder,result_file,threshold_in_point):
    def sum_sta(day_folder):
        sta_list=[sac_file[:-4] for sac_file in os.listdir(day_folder)]
        return list(set(sta_list))

    #print (day_folder)
    sta_list = sum_sta(os.path.join(year_folder, day_folder))

    delete_sta_list = []
    delete_sta_list_append = delete_sta_list.append
    for sta in sta_list:
        files = glob.glob(os.path.join(year_folder, day_folder, ''.join(['*', sta, '*'])))
        if len(files) < 3:
            delete_sta_list_append(sta)
            continue
        st = obspy.read(os.path.join(year_folder, day_folder, ''.join(['*', sta, '*'])))
        for tr in st:
            if checkInterrupt(tr, threshold_in_point=threshold_in_point):
                delete_sta_list_append(sta)
                break

    if delete_sta_list != []:
        with open(result_file, 'a') as f:
            f.write(' '.join([str(year_folder), str(day_folder), str(delete_sta_list), '\n']))
def checkLack_Interrupt_recordOneday_parallel(year_folder,result_file,threshold_in_point=6000):

    cores = multiprocessing.cpu_count()
    pool_new = multiprocessing.Pool(processes=5)

    tasks=[]
    for day_folder in sorted(os.listdir(year_folder)):
        tasks.append((day_folder,year_folder,result_file,threshold_in_point))
    #print (len(tasks))
    rs = pool_new.starmap_async(checkInterrupt_day, tasks,chunksize=1)  # chunksize is how many tasks will be processed by one processor
    # close() & join() is necessary
    pool_new.close()
    # simple progress bar
    while (True):
        if (rs.ready()): break
        remaining = rs._number_left
        print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
              end='\r')  # '\r' means remove the last line
        time.sleep(0.5)
    #print ('pool')
    pool_new.join()

def checkLack_recordOneday(year_folder,result_file):
    def sum_sta(day_folder):
        sta_list=[]
        for sac_file in os.listdir(day_folder):
            sta=sac_file[:-4]
            if sta not in sta_list:
                sta_list.append(sta)
        return sta_list


    for day_folder in sorted(os.listdir(year_folder)):
        sta_list=sum_sta(os.path.join(year_folder,day_folder))

        delete_sta_list = []
        for sta in sta_list:
            lack=False
            files = glob.glob(os.path.join(year_folder, day_folder,'*'+sta+'*'))
            if len(files) < 3:
                lack=True

            if lack:
                delete_sta_list.append(sta)
                with open(result_file, 'a') as f:
                    f.write(str(year_folder) + ' ' + str(day_folder) + ' ' + sta + '\n')
        logging.info('delete the data of these stations on %s: %s'%(day_folder,str(delete_sta_list)))
if __name__ == '__main__':
    file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/RAW_data/allData_test/2010/20100119/NC.GDXB.HHZ'
    st=obspy.read(file)
    trace=st[0]

    start = time.time()
    result=checkInterrupt(trace)
    end=time.time()
    print (end-start)
    print (result)
    print ('!!!!!!!')

