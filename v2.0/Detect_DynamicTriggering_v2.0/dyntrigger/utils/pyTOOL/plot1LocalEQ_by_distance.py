"""
Written by Naidan YUN on 20180523.
"""
import os
import obspy
import matplotlib.pyplot as plt

DATA_FOLDER='/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang/event_data/removeRES/events'
EVENT='20130208152638.470'

def plot_run(tb,te,t_name):
    st=obspy.read(os.path.join(DATA_FOLDER,EVENT,'*BHZ'))
    st.detrend('linear')
    st.detrend('constant')
    st.filter('highpass', freq=5)
    delta=st[0].stats.delta
    b=st[0].stats.sac['b']
    indexB = int((tb-b)/delta)
    indexE = int((te-b)/delta)
    p_data={}
    p_sta={}
    for tr in st:
        if t_name in tr.stats.sac:
            p_data[tr.stats.sac[t_name]]=tr.data[indexB:indexE]
            p_sta[tr.stats.sac[t_name]] = tr.stats.station
    p_data=sorted(p_data.items(),key=lambda item:item[0])
    p_sta = sorted(p_sta.items(), key=lambda item: item[0])
    fig=plt.figure(figsize=[8,12])
    #ax1=fig.add_subplot(111)
    for i in range(len(p_data)):
        ax1 = fig.add_subplot(len(p_data),1,i+1)
        tar_data=p_data[i][1]
        sta_name=p_sta[i][1]
        #ax1.plot([i for i in range(len(tar_data))],tar_data/2+i*2000,'black')
        p_time=p_data[i][0]
        p_index=int((p_time-b)/delta)
        ax1.plot([i for i in range(len(tar_data))], tar_data, 'black')
        ax1.plot([p_index-indexB,p_index-indexB],[min(tar_data),max(tar_data)],'red')
        ax1.set_yticks([0])
        ax1.set_yticklabels([sta_name])
        ax1.set_xticks([1000 * i for i in range(((indexE - indexB) // 1000)+1)])
        if i !=0:
            ax1.spines['top'].set_color('none')
        if i != len(p_data)-1:
            ax1.set_xticks([])
            ax1.spines['bottom'].set_color('none')
        if i == len(p_data)-1:
            ax1.set_xticklabels([int((indexB + 1000 * i)/100+b) for i in range(((indexE - indexB) // 1000)+1)])
            ax1.set_xlabel('Time (s)')

    return fig

if __name__ == "__main__":
    fig=plot_run(2480,2530,'t0')
    fig.show()
    fig.savefig('/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Result/Xiaojiang/Belta_test/Spec/20130208/manul_pick.pdf')
    print ('Finish')





