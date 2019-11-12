"""
Written by Naidan YUN on 20180721.
"""
import os
import subprocess
import time
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')

datapath= '/data4/yunnd/rawData_100HZ/2013'
g = os.walk(datapath)

PZfile='/data1/yunnd/Dynamic_triggering/Xiaojiang/data/removeResponse/PZfile/QJPZ_single.txt'

outpath='/data4/yunnd/removeRES_100HZ/2013'

for path,dir_list,file_list in g:
    for file_name in sorted(file_list):
        start=time.time()
        input_file=os.path.join(path, file_name)

        output_path=os.path.join(outpath,path.replace(datapath+'/',''))
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_file=os.path.join(output_path,file_name)
        print (input_file)
        print (output_path)
        p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
        s = "wild echo off \n"
        s += "r %s \n"%input_file
        s += "ch LCALDA TRUE \n"
        s += "ch LOVROK TRUE \n"
        s += "rmean;rtr;taper \n"
        s += "trans from polezero subtype %s to vel freq 0.005 0.01 30 35\n"%PZfile
        s += "mul 1.0e9 \n"
        s += "rmean;rtr \n"
        s += "w %s \n"%output_file
        s += "q \n"
        p.communicate(s.encode())
        end=time.time()
        print ('Using time: %f'%(end-start))
print ('fnish')
