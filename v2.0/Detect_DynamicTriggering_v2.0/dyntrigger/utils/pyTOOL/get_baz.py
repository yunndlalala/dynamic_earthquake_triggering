import os
import obspy
import pandas as pd
out_file='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/tele_baz.csv'
data_folder='C:/Users/yunnd/Documents/Research/Dynamic_triggering/Xiaojiang/data/events_v2/events'
tele_baz=[]
for event in os.listdir(data_folder):
    print (event)
    event_path=os.path.join(data_folder,event)
    first_sac_file=os.listdir(event_path)[0]
    sac_file_path=os.path.join(event_path,first_sac_file)
    st=obspy.read(sac_file_path)
    tr=st[0]
    baz=tr.stats.sac['baz']
    tele='%s-%s-%sT%s:%s:%sZ'%(event[:4],event[4:6],event[6:8],event[8:10],event[10:12],event[12:])
    tele_baz.append([tele,baz])

out_df=pd.DataFrame(data=tele_baz,columns=['event','baz'])
out_df.to_csv(out_file,index=False)
print ('Finish!')