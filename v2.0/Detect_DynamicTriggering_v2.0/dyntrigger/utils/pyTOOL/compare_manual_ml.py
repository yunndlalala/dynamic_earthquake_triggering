'''
Written by Naidan YUN on 20180323.
This script is used to compare the hand-picked eqs and the catalog
which is the result of Match&Locte in order to check whether the latter contains the former.
And the catalog of hand-picked eqs is from the script 'get_phase_time.py'.
'''
import os
import glob
import re
import numpy as np
import pandas as pd
from datetime import datetime,timedelta
from TOOL import datetime2UTCstr

manual_eq_file='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/local_phase_time.csv'
result_ML_path='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/Match_Locate/Result_of_lyx/DetectedFinal'
out_file='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/result_compare_manual_ML.csv'

def preprocess_reult_ML(result_ML_path):
    #preprocessing the result of ML in order to get a whole catalog
    day_files=glob.glob(os.path.join(result_ML_path,'*/*'))
    catalog=[]
    for day_file in day_files:
        with open(day_file) as f:
            data=f.readlines()[1:]
            original_time=[datetime.strptime(event[7:32],'%Y/%m/%d %H:%M:%S.%f') for event in data]
            catalog=catalog+original_time
    return catalog

def preprocess_result_manual(manual_eq_file):
    #preprocessing the result of manual eqs in order to get a whole catalog
    manual_eq = pd.read_csv(manual_eq_file)
    manual_eq['ot'] = pd.to_datetime(manual_eq['ot'])
    catalog=manual_eq.groupby(['teleseism','local_eq_num'])['ot'].min()
    return catalog

if __name__ == "__main__":
    #compare the above two catalogs
    manual_catalog=preprocess_result_manual(manual_eq_file) #dataframe
    ML_catalog=preprocess_reult_ML(result_ML_path) #list
    min_dt=[]
    mark=[]
    corresponding_ML_eq=[]
    for manual_eq_ot in manual_catalog:
        print ('The manual eq is :',manual_eq_ot)
        flag = False
        #tar_ML_eq_ot='Null'
        abs_seconds=[]
        for ML_eq_ot in ML_catalog:
            abs_seconds.append(abs((ML_eq_ot - manual_eq_ot).total_seconds()))
        tar_ML_eq_ot=ML_catalog[abs_seconds.index(min(abs_seconds))]
        tar_ML_eq_ot = datetime2UTCstr(tar_ML_eq_ot)
        print('The corresponding eq in ML results is : ', tar_ML_eq_ot)
        if min(abs_seconds)< 600:
            flag=True
        print (flag)
        min_dt.append(min(abs_seconds))
        mark.append(flag)
        corresponding_ML_eq.append(tar_ML_eq_ot)
    out_frame=manual_catalog.to_frame()
    out_frame['ot']=[datetime2UTCstr(i) for i in out_frame['ot']]
    out_frame['min_dt'] = min_dt
    out_frame['mark']=mark
    out_frame['corresponding_ML_eq']=corresponding_ML_eq
    out_frame.to_csv(out_file)


    print ('Finish')

















