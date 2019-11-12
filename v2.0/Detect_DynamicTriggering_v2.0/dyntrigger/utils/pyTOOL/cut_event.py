"""
cutevent by Yijian ZHOU 2018-01-02
modified by Naidan YUN 2018-03-13
input:
catalog; stream_path; out_path; time window
output:
./events/[origin time]/net.sta.ot.chn.SAC
head file:
stlo stla stel
evlo evla evdp
kztime
b = 0
other:
The structure of data storage must be [year]/[yearmonthday]/sac file
The catalog must be
                    time,latitude,longitude,depth,mag
                    Y-m-dTH:M:SZ,lat,lon,depth,mag
                    .............................
                    .............................


Modified by Naidan YUN on 20180706.
    Fixed a problem with the time variable in the header.
    Turn the "print" to "logging" which is more normalize.
"""
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
import os, glob, shutil
from obspy.core import *
import subprocess
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')

import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--stream_path", type = str,
                    default = '/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Geysers/Data/removeRES/GDX/all_data')
parser.add_argument("--catalog", type = str,
                    default = '/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Geysers/Data/catalog4cut_event_Hifi_tri.csv')
parser.add_argument("--same_time_win", type = str,
                    default='False')
parser.add_argument("--time_win_before", type = float,
                    default=-18000)
parser.add_argument("--time_win_after", type = float,
                    default=18000)
parser.add_argument("--out_path", type = str,
                    default ='/Users/yunnaidan/Project/Dynamic_Triggering/Workspace/Geysers/Data/removeRES/GDX/hifi_tri_events')
args = parser.parse_args()

# set the time window
same_time_win=args.same_time_win
time_before = args.time_win_before
time_after = args.time_win_after
# i/o path
stream_path = args.stream_path
out_path = os.path.join(args.out_path,'events')
if os.path.exists(out_path):
    raise ValueError('check your output path!')
else:
    os.mkdir(out_path)

tmp = open(args.catalog); ctlg_lines = tmp.readlines(); tmp.close()
for ctlg_line in ctlg_lines[1:]:
    print('cutting event {}'.format(ctlg_line))
    #print (ctlg_line)
    if same_time_win=='True':
        date_time, lat, lon, depth, mag = ctlg_line.split(',')
    else:
        date_time, lat, lon, depth, mag,time_before_single, time_after_single= ctlg_line.split(',')
        time_before_single=float(time_before_single)
        time_after_single=float(time_after_single)
        time_before=time_before_single
        time_after=time_after_single

    mag = mag.split('\n')[0]
    t0 = UTCDateTime(date_time)
    # make output dir
    date, time = date_time.split('T')
    time_key = ''.join(date.split('-')) +\
               ''.join(time.split(':'))[:-1]

    # time window for slicing

    ts = t0 + time_before
    te = t0 + time_after

    # find stream paths
    rela_path = ''.join([str(ts.year),'%02d'%ts.month,'%02d'%ts.day])
    path = os.path.join(stream_path, str(ts.year), rela_path)
    # use sac
    if os.path.exists(path):
        os.chdir(path)
        out_dir = os.path.join(out_path, time_key)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for stream in glob.glob('*'):
            print (stream)
            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
            s = "wild echo off \n"
            #print (stream)
            net, sta, chn= stream.split('.')
            fname = '.'.join([net, sta, chn])
            #print (fname)
            b = ts - UTCDateTime('%04d'%ts.year+'%02d'%ts.month+'%02d'%ts.day)
            e = te - UTCDateTime('%04d'%ts.year+'%02d'%ts.month+'%02d'%ts.day)
            # cut event and change the head file
            flag=False
            if te.day != ts.day:
                next_stream=os.path.join(stream_path,str(te.year),
                                         ''.join([str(te.year),'%02d'%te.month,'%02d'%te.day]),stream)
                if os.path.exists(next_stream):
                    flag=True
                    s += "r %s \n" % os.path.join(path, stream)
                    s += "merge GAP ZERO o a %s \n"%next_stream
                    #print (os.path.join(path, 'long.'+stream))
                    s += "w %s \n" % (os.path.join(path, 'long.'+stream))
                    stream='long.'+stream
                else:
                    logging.warning('The data of next day is not found!')
            #print (flag)
            #print('cutting stream {}, {} to {}'.format(stream, b, e))
            s += "cuterr fillz \n"
            s += "cut b %s %s \n" %(b, e)
            s += "r %s \n" %os.path.join(path, stream)
            s += "ch LCALDA TRUE \n"
            s += "ch LOVROK TRUE \n"
            #s += "ch IZTYPE IO  \n"
            s += "ch b %s \n"%str(time_before)
            s += "ch e %s \n" % str(time_after)
            s += "ch nzyear %s nzjday %s \n" %(str(t0.year),str(t0.julday))
            s += "ch nzhour %s nzmin %s nzsec %s \n" %(str(t0.hour), str(t0.minute), str(t0.second))
            s += "ch nzmsec %s \n" %str(t0.microsecond)[0:3]
            s += "ch evlo %s evla %s evdp %s \n" %(lon, lat,depth)
            s += "ch mag %s \n" %mag
            #print ('!!!!!!!!!',fname)
            s += "w %s \n" %os.path.join(out_dir, fname)
            #print (flag)
            s += "q \n"
            p.communicate(s.encode())
            if flag:
                os.remove(os.path.join(path, stream))
    else:
        logging.warning('%s is not exist!'% path)


