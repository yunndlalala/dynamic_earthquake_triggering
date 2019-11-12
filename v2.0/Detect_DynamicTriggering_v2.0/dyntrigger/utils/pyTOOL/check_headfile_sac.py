"""
Written by Naidan YUN on 20180821.
In order to check whether the wave data has been removed response.

Importantly, the parallel doesn't have as much improvement as other script for running, maybe because the task is I/O intensive.
And using the order 'saclst' is very much faster than using obspy.
"""
import obspy
import os
import multiprocessing
import time
import re
def check_1file_IDEP(file,out_file):
    st = obspy.read(file)
    tr = st[0]
    if 'idep' not in tr.stats['sac'].keys():
        with open(out_file, 'a') as f:
            f.write(file+'\n')
    else:
        IDEP = tr.stats['sac']['idep']
        if IDEP != 7:
            with open(out_file, 'a') as f:
                f.write(file + '\n')
def check_1file_IDEP_sac(file,out_file):
    result = os.popen('saclst IDEP f %s'%file)
    for line in result:
        value = re.split(' +', line)[1][0]
    if value == '-12345':
        with open(out_file, 'a') as f:
            f.write(file + '\n')
    else:
        if value != '7':
            with open(out_file, 'a') as f:
                f.write(file + '\n')
def check_IDEP(data_path,out_file):
    if os.path.exists(out_file):
        os.remove(out_file)
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    tasks=[]
    g=os.walk(data_path)
    for path,dir_list,file_list in g:
        for file in sorted(file_list):
            tasks.append((os.path.join(path,file),out_file))
    rs=pool.starmap_async(check_1file_IDEP_sac, tasks,chunksize=1)
    pool.close()

    while (True):
        if (rs.ready()): break
        remaining = rs._number_left
        print ("finished:{0}/{1}".format(len(tasks)-remaining,len(tasks)),end='\r')#'\r' means remove the last line
        time.sleep(0.5)

    pool.join()
if __name__ == '__main__':
    data_path = '/data4/yunnd/removeRES/For_cal_backgroundNoisePID'
    #data_path='/data4/yunnd/RAW/rawData_100HZ_aboutTeleseism/2012/20120319'
    outfile='checkIDEPResult_of_For_cal_backgroundNoisePID.txt'
    check_IDEP(data_path,out_file=outfile)
    print ('Finish!!!!!!!!!!!!!!!!')

