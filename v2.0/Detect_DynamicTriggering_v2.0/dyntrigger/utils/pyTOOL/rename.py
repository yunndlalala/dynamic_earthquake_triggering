"""
@version:
author:yunnaidan
@time: 2019/07/07
@file: rename.py
@function:
"""
import os
data_path='/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Dynamic_triggering_project/Kunming/data/raw/2018'
for day_folder in os.listdir(data_path):
    if day_folder != '.DS_Store':
        for sacfile in os.listdir(os.path.join(data_path,day_folder)):
            if sacfile !='.DS_Store':
                print (sacfile)
                net=sacfile.split('.')[6]
                sta = sacfile.split('.')[7]
                chn=sacfile.split('.')[9]
                new_name='.'.join([net,sta,chn])
                os.system('mv %s %s'%(os.path.join(data_path,day_folder,sacfile),os.path.join(data_path,day_folder,new_name)))
print ('finish')
