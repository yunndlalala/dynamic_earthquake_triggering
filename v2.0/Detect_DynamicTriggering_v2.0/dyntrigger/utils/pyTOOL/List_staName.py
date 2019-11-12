"""
Written by Naidan YUN on 20180722.
"""
import os
import shutil


def list_staName(file_dir, outfile):
    sta_name_list = []
    g = os.walk(file_dir)
    for path, dir_list, file_list in g:
        for file in file_list:
            sta_name = file.split('.')[1]
            if sta_name not in sta_name_list:
                sta_name_list.append(sta_name)
    with open(outfile, 'w') as f:
        for sta_name in sorted(sta_name_list):
            f.write(sta_name + '\n')


if __name__ == "__main__":
    file_dir = 'C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/events_v2/events'
    outfile = 'C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/events_v2/contain_stations.txt'
    list_staName(file_dir, outfile)

    print('Finish')
