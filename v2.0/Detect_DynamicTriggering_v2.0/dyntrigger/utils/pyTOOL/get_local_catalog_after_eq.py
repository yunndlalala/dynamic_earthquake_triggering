'''
20180318 Naidan YUN
This script is used to select the local earthquakes which happened
in the 9000s after the distant earthquake. So a local earthquake catalog
is needed in this process.
As far, the catalog is from Chuan YAN, which is used to be the template in Match&Locate.
'''
import pandas as pd
from datetime import datetime,timedelta

DISTANT_EQ_CATALOG='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/EQ_info.csv'
LOCAL_EQ_CATALOG='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/Match_Locate/catalog.dat'
TIME_WINDOW_LENGTH=9000
OUT_FILE='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/local_earthquakes.dat'

#load the catalog of distant eq
dis_cata_data=pd.read_csv(DISTANT_EQ_CATALOG,encoding='gb2312')
dis_cata=[datetime.strptime(i[:10]+' '+i[11:-1],'%Y-%m-%d %H:%M:%S.%f') for i in dis_cata_data['time'].values]
#load the catalog of local eq
with open(LOCAL_EQ_CATALOG) as f:
    loc_cata_data=f.read()
    loc_cata_data=loc_cata_data.split('\n')
    loc_cata=[datetime.strptime(i.split(' ')[0]+' '+i.split(' ')[1],'%Y/%m/%d %H:%M:%S.%f') for i in loc_cata_data]

#begin search
out_eqs={}
for i in range(len(dis_cata)):
    eq_time=datetime.strftime(dis_cata[i],'%Y-%m-%d'+'T'+'%H:%M:%S.%f'+'Z')
    print (eq_time)
    out_eqs[eq_time]=[]
    for j in range(len(loc_cata)):
        if loc_cata[j] >=dis_cata[i] and  loc_cata[j] <= dis_cata[i]+timedelta(seconds=TIME_WINDOW_LENGTH):
            out_eqs[eq_time].append(loc_cata_data[j])

#write the result
with open(OUT_FILE,'w') as out_file:
    out_file.write(str(out_eqs))


print ('Finish!!!!!!!!')




