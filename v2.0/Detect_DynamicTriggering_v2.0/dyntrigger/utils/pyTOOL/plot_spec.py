"""
@version:
author:yunnaidan
@time: 2018/10/21
@file: plot_spec.py
@function:
"""
import numpy as np
import matplotlib.pyplot as plt

import obspy
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt
from obspy import UTCDateTime

def spec(st,ax,begin_t,end_t):

    tr = st[0]
    npts = tr.stats.npts
    dt = tr.stats.delta
    t = np.linspace(begin_t,end_t, (end_t-begin_t)/dt)
    f_min = 1
    f_max = 50

    begin_point=begin_t/dt
    end_point=end_t/dt
    scalogram = cwt(tr.data[begin_point:end_point], dt, 8, f_min, f_max)


    x, y = np.meshgrid(
        t,
        np.logspace(np.log10(f_min), np.log10(f_max), scalogram.shape[0]))

    ax.pcolormesh(x, y, np.abs(scalogram), cmap=obspy_sequential)
    ax.set_xlabel("Time after %s [s]" % tr.stats.starttime)
    ax.set_ylabel("Frequency [Hz]")
    ax.set_yscale('log')
    ax.set_ylim(f_min, f_max)
    return ax
if __name__=='__main__':
    sacfile='/Users/yunnaidan/XJ.QJ01.BHE'
    st = obspy.read(sacfile)
    ot = UTCDateTime("2012-03-20T14:18:54")
    t1=ot+14400
    t2=t1+9000
    st=st.slice(t1,t2)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    st.spectrogram(log=True,axes=ax1,show=False)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    spec(st, ax)
    fig.show()
