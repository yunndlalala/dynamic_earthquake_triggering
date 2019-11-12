"""
Written by Naidan YUN on 20180722.
"""
import os
import shutil
import obspy
import numpy as np
def copy_file_and_dir(datapath,sta_list,outpath):
    g = os.walk(datapath)
    for path, dir_list, file_list in g:
        for file in file_list:
            sta_name = file.split('.')[1]
            if sta_name in sta_list:
                print (file)
                outdir=os.path.join(outpath,path.replace(datapath+'/',''))
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                shutil.copy(os.path.join(path,file),outdir)
def compare_sacfile_data(sacfile1,sacfile2):
    tr1=obspy.read(sacfile1)[0]
    tr2=obspy.read(sacfile2)[0]
    data1=tr1.data
    data2=tr2.data
    compare_result=data1-data2
    compare_mean=np.mean(compare_result)
    return compare_mean

if __name__ == "__main__":
    # datapath=''
    # sta_list=[]
    # outpath=''
    # copy_file_and_dir(datapath, sta_list, outpath)
    sacfile1='/data4/yunnd/removeRES/removeRESData_100HZ_aboutTeleseism/2012/20120319/XJ.QJ01.BHN'
    sacfile2='/data4/yunnd/removeRES/test_remove_correctness/20120319/XJ.QJ01.BHN'
    print (compare_sacfile_data(sacfile1,sacfile2))
    print ('finish')