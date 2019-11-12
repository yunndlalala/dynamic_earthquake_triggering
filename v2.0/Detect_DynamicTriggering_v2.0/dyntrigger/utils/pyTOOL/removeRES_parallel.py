"""
Written by Naidan YUN on 20180721.
"""
import os
import subprocess
import time
import multiprocessing
from multiprocessing import Lock
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')



def removeRES_for_1file_throughPZfile_parallel(file_name,input_path,output_path,PZfile):
    lock.acquire()
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    lock.release()
    output_file = os.path.join(output_path, file_name)
    input_file= os.path.join(input_path, file_name)
    print(input_file)
    print (output_file)
    #print(output_file)
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "r %s \n" % input_file
    s += "ch LCALDA TRUE \n"
    s += "ch LOVROK TRUE \n"
    #s += "rmean;rtr;taper \n"
    s += "rmean;rtr\n"
    s += "trans from polezero subtype %s to vel freq 0.005 0.01 40 45\n" % PZfile # for broadband
    #s += "trans from polezero subtype %s to vel freq 0.1 0.5 40 45\n" % PZfile  # for short period station
    s += "mul 1.0e9 \n"
    s += "rmean;rtr \n"
    s += "w %s \n" % output_file
    s += "q \n"
    p.communicate(s.encode())
# def removeRES_for_1file_throughRESPfile_parallel(file_name, input_path, output_path, RESPfile):
#     lock.acquire()
#     if not os.path.exists(output_path):
#         os.makedirs(output_path)
#     lock.release()
#     output_file = os.path.join(output_path, file_name)
#     input_file= os.path.join(input_path, file_name)
#     print(input_file)
#     #print(output_file)
#     p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
#     s = "wild echo off \n"
#     s += "r %s \n" % input_file
#     s += "ch LCALDA TRUE \n"
#     s += "ch LOVROK TRUE \n"
#     s += "rmean;rtr;taper \n"
#     s += " trans from evalresp fname %s to vel freq 0.005 0.01 30 35\n" % RESPfile
#     #s += "mul 1.0e9 \n"
#     s += "rmean;rtr \n"
#     s += "w %s \n" % output_file
#     s += "q \n"
#     p.communicate(s.encode())
#
# def removeRES_for_1file_throughRESPfile_notaper(file_name, input_path, output_path, RESPfile):
#     if not os.path.exists(output_path):
#         os.makedirs(output_path)
#     output_file = os.path.join(output_path, file_name)
#     input_file= os.path.join(input_path, file_name)
#     print(input_file)
#     #print(output_file)
#     p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
#     s = "wild echo off \n"
#     s += "r %s \n" % input_file
#     s += "ch LCALDA TRUE \n"
#     s += "ch LOVROK TRUE \n"
#     s += "rmean;rtr \n"
#     s += " trans from evalresp fname %s to vel freq 0.005 0.01 30 35\n" % RESPfile
#     #s += "mul 1.0e9 \n"
#     s += "rmean;rtr \n"
#     s += "w %s \n" % output_file
#     s += "q \n"
#     p.communicate(s.encode())
def removeRES_for_1file_throughPZfile(file_name,input_path,output_path,PZfile):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file = os.path.join(output_path, file_name)
    print (output_file)
    input_file= os.path.join(input_path, file_name)
    print(input_file)
    #print(output_file)
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "r %s \n" % input_file
    s += "ch LCALDA TRUE \n"
    s += "ch LOVROK TRUE \n"
    #s += "rmean;rtr;taper \n"
    s += "rmean;rtr\n"
    s += "trans from polezero subtype %s to vel freq 0.005 0.01 30 35\n" % PZfile
    s += "mul 1.0e9 \n"
    s += "rmean;rtr \n"
    s += "w %s \n" % output_file
    s += "q \n"
    p.communicate(s.encode())
def removeRES_for_1file_throughPZfile_simple(input_file,output_file,PZfile):
    #print(output_file)
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "r %s \n" % input_file
    s += "ch LCALDA TRUE \n"
    s += "ch LOVROK TRUE \n"
    s += "rmean;rtr \n"
    s += "trans from polezero subtype %s to vel freq 0.005 0.01 30 35\n" % PZfile
    s += "mul 1.0e9 \n"
    s += "rmean;rtr \n"
    s += "w %s \n" % output_file
    s += "q \n"
    p.communicate(s.encode())
def init(l):
    global lock
    lock = l

def remove_res(datapath, resfiles_dict, outdir):
    lock = Lock()
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=15,initializer=init, initargs=(lock,))

    g = os.walk(datapath)

    tasks=[]
    for path,dir_list,file_list in g:
        for file in file_list:
            sta='.'.join(file.split('.')[:2])
            chn=file.split('.')[2]
            PZfile=resfiles_dict[(sta,chn)]

            tasks.append((file, path, os.path.join(outdir, path.split('/')[-1]), PZfile))


    rs=pool.starmap_async(removeRES_for_1file_throughPZfile_parallel, tasks, chunksize=1)

    pool.close()

    while (True):
        if (rs.ready()): break
        remaining = rs._number_left
        print ("finished:{0}/{1}".format(len(tasks)-remaining,len(tasks)),end='\r')#'\r' means remove the last line
        time.sleep(0.5)

    pool.join()
# def main_Xiaojiang():
#     lock = Lock()
#     cores = multiprocessing.cpu_count()
#     pool = multiprocessing.Pool(processes=14,initializer=init, initargs=(lock,))
#
#     datapath = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/event_data/RAW_data/events'
#     g = os.walk(datapath)
#     PZfile = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/removeResponse/QJPZ_single.txt'
#     outdir = '/data1/yunnd/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/event_data/removeRES/events'
#
#     tasks=[]
#     for path,dir_list,file_list in g:
#         for file in sorted(file_list):
#             tasks.append((file,path,os.path.join(outdir,path.replace(datapath+'/','')),PZfile))
#
#
#     rs=pool.starmap_async(removeRES_for_1file_throughPZfile_parallel, tasks, chunksize=1)
#
#     pool.close()
#
#     while (True):
#         if (rs.ready()): break
#         remaining = rs._number_left
#         print ("finished:{0}/{1}".format(len(tasks)-remaining,len(tasks)),end='\r')#'\r' means remove the last line
#         time.sleep(0.5)
#
#     pool.join()
if __name__ =="__main__":
    #main_Xiaojiang()
    print ('finish!!!!!!!!!!!!!')
