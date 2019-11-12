"""
@version:
author:yunnaidan
@time: 2018/09/10
@file: compare_single_parallel.py
@function:
"""
import numpy as np
single_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/result/summary_14400_removeRES.txt'
parallel_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/result/summary_14400_removeRES_Singleprocess.txt'

with open(single_file) as f:
    data1=f.readlines()
data1_dict={}
for line in data1:
    tele=line.split(' ')[0]
    sta=line.split(' ')[1]
    key=(tele,sta)
    value=[line.split(' ')[2],line.split(' ')[3][:-1]]
    data1_dict[key]=value

with open(parallel_file) as f:
    data2=f.readlines()
data2_dict={}
for line in data2:
    tele=line.split(' ')[0]
    sta=line.split(' ')[1]
    key=(tele,sta)
    value=[line.split(' ')[2],line.split(' ')[3][:-1]]
    data2_dict[key]=value

com=[]
for key in data1_dict.keys():
    value1=data1_dict[key]
    value2=data2_dict[key]
    com0=abs(float(value1[0])-float(value2[0]))
    com1=abs(float(value1[1])-float(value2[1]))
    com.append(com0)
    com.append(com1)

print (np.mean(com))
print ('finish')