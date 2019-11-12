"""
Written by Naidan YUN on 20180715.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Info_file='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/data/Event_allResult_Info.csv'
dataframe=pd.read_csv(Info_file)
DS_Z=dataframe['DS_Z'].values*(10**(-6))
DS_H=dataframe['DS_H'].values*(10**(-6))
HighPID=dataframe['HighPID'].values
belta_trigger=dataframe['belta_trigger'].values
HighPID_trigger=dataframe['HighPID_trigger'].values
HighPID_log=[]
for item in HighPID:
    if item > 0:
        HighPID_log.append(np.log10(item))
    else:
        HighPID_log.append(-np.log10(-item))
fig=plt.figure()
ax1=fig.add_subplot(211)
for i in range(len(DS_Z)):
    if belta_trigger[i]==1 and HighPID_trigger[i]==1:
        ax1.scatter(DS_Z[i], HighPID_log[i], s=120, color='r', edgecolors='')
    elif belta_trigger[i] ==1 and HighPID_trigger[i] == 0:
        ax1.scatter(DS_Z[i], HighPID_log[i], s=120, color='orange', edgecolors='')
    elif belta_trigger[i] ==0 and HighPID_trigger[i] == 1:
        ax1.scatter(DS_Z[i], HighPID_log[i], s=120, color='green', edgecolors='')
    elif belta_trigger[i] ==0 and HighPID_trigger[i] == 0:
        ax1.scatter(DS_Z[i], HighPID_log[i], s=120, color='lightgray', edgecolors='')

ax1.set_xscale('log')
ax1.set_yticks([-6,-4,-2,0,2,4,6,8])
ax1.set_yticklabels(['-10^6','-10^4','-10^2','0','10^2','10^4','10^6','10^8'])
ax1.set_xlim([4*10**(-5),6*10**(-2)])
plt.xlabel('DS (Mpa)')
plt.ylabel('Power integral/(nm/s)^2')

ax2=fig.add_subplot(212)
for i in range(len(DS_Z)):
    if belta_trigger[i]==1 and HighPID_trigger[i]==1:
        ax2.scatter(DS_H[i], HighPID_log[i], s=120, color='r', edgecolors='')
    elif belta_trigger[i] ==1 and HighPID_trigger[i] == 0:
        ax2.scatter(DS_H[i], HighPID_log[i], s=120, color='orange', edgecolors='')
    elif belta_trigger[i] ==0 and HighPID_trigger[i] == 1:
        ax2.scatter(DS_H[i], HighPID_log[i], s=120, color='green', edgecolors='')
    elif belta_trigger[i] ==0 and HighPID_trigger[i] == 0:
        ax2.scatter(DS_H[i], HighPID_log[i], s=120, color='lightgray', edgecolors='')
ax2.set_xscale('log')
ax2.set_yticks([-6,-4,-2,0,2,4,6,8])
ax2.set_yticklabels(['-10^6','-10^4','-10^2','0','10^2','10^4','10^6','10^8'])
ax2.set_xlim([4*10**(-5),6*10**(-2)])
plt.xlabel('DS (Mpa)')
plt.ylabel('Power integral/(nm/s)^2')
fig.show()

print ('finish!')