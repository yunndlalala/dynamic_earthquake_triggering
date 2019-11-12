"""
@version:
author:yunnaidan
@time: 2018/09/12
@file: remove_QJ19_from_summaryfile.py
@function:
"""
out_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/result/summary_14400_removeRES_mag6_7.txt'
filter_out_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Energy_Ratio/result/summary_14400_removeRES_mag6_7_filter.txt'
with open(out_file) as f_out_file:
    already_data = f_out_file.readlines()
filter_data = []
for line in already_data:
    sta=line.split(' ')[1]
    b_PI=float(line.split(' ')[2])
    a_PI = float(line.split(' ')[3][:-1])
    if b_PI > 1 and a_PI > 1 and sta!='QJ19':
        filter_data.append(line)

with open(filter_out_file,'a') as f_filter_out_file:
    for line in filter_data:
        f_filter_out_file.write(line.split(' ')[0] + ' ' +
                                line.split(' ')[1] + ' ' +
                                line.split(' ')[2] + ' ' +
                                line.split(' ')[3][:-1]+'\n')

print ('finish')

