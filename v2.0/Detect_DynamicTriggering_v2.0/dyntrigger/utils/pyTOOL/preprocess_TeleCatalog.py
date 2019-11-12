import pandas as pd
import numpy as np
import os
import obspy
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from obspy.taup import TauPyModel
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
model = TauPyModel(model="iasp91")
from dyntrigger.utils.pyTOOL.cal_distance import haversine

# TOOL


def str2datetime(str_list):
    datetime_lst = [datetime.strptime(
        i[:10] + ' ' + i[11:-1], '%Y-%m-%d %H:%M:%S.%f') for i in str_list]
    return datetime_lst


def a_abstime(ot, depth, distance_in_degree):
    arrivals = model.get_travel_times(source_depth_in_km=depth,
                                      distance_in_degree=distance_in_degree)
    B_time = ot + arrivals[0].time
    return B_time


def vel_abstime(ot, vel, distance_in_km):
    travel_time = distance_in_km / vel
    E_time = ot + travel_time
    return E_time


def make_stalta_as_Etime_abs(
        trace,
        ot,
        sta,
        lta,
        thr_on,
        thr_off,
        fig_outfile=None):
    df = trace.stats.sampling_rate
    tar_trace = trace.slice(
        starttime=ot,
        endtime=ot + 7200,
        nearest_sample=True)
    cft = classic_sta_lta(tar_trace, int(sta * df), int(lta * df))
    on_off = np.array(trigger_onset(cft, thr_on, thr_off))
    if len(on_off) == 0:
        Etime = 0
    else:
        Epoint = on_off[0][1]
        Etime = ot + Epoint / df

    if fig_outfile is not None:
        fig = plot_trigger(tar_trace, cft, thr_on, thr_off, show=False)
        fig.savefig(fig_outfile)
    return Etime

# Main functions for time windows


def BE_A_vel(dat_path, dat_file):
    dat = pd.read_csv(os.path.join(dat_path, dat_file))
    # calculate B_time and E_time
    B_time_list = []
    E_time_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        dist_km = dat.iloc[i]['dist(km)']
        distance_in_degree = dist_km / 111
        B_time = a_abstime(ot, depth, distance_in_degree)

        vel = 2
        E_time = vel_abstime(ot, vel, dist_km)

        B_time_list.append(str(B_time))
        E_time_list.append(str(E_time))
    dat['B_time'] = B_time_list
    dat['E_time'] = E_time_list

    dat.to_csv(os.path.join(dat_path, dat_file), index=False)
    return None


def Bvel_Evel(catalog_file, sta_lat,sta_lon,Tb, vel_B, vel_E):
    dat = pd.read_csv(catalog_file)
    # calculate B_time and E_time
    A_time_list = []
    Tb_B_list = []
    Tb_E_list = []
    Te_B_list = []
    Te_E_list = []
    dist_list=[]
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        event_lat = dat.iloc[i]['latitude']
        event_lon = dat.iloc[i]['longitude']
        if sta_lon == None and sta_lat == None:
            dist_km = dat.iloc[i]['dist(km)']
            distance_in_degree = dist_km / 111
        else:
            dist_km,distance_in_degree=haversine(sta_lon,sta_lat,event_lon,event_lat)


        A_time = a_abstime(ot, depth, distance_in_degree)
        Tb_B = A_time - Tb
        Tb_E = A_time
        Te_B = vel_abstime(ot, vel_B, dist_km)
        Te_E = vel_abstime(ot, vel_E, dist_km)


        A_time_list.append(str(A_time))
        Tb_B_list.append(str(Tb_B))
        Tb_E_list.append(str(Tb_E))
        Te_B_list.append(str(Te_B))
        Te_E_list.append(str(Te_E))
        dist_list.append(dist_km)

    dat['A_time'] = A_time_list
    dat['Tb_B_time'] = Tb_B_list
    dat['Tb_E_time'] = Tb_E_list
    dat['Te_B_time'] = Te_B_list
    dat['Te_E_time'] = Te_E_list
    if sta_lon != None and sta_lat != None:
        dat['dist(km)'] = dist_list

    dat.to_csv(catalog_file, index=False)

    return None


def BE_stalta(
        catalog_file,
        data_path,
        sta_lta_parameter,
        catalog_outfile,
        fig_outpath=None):
    dat = pd.read_csv(catalog_file)
    # calculate B_time and E_time
    B_time_list = []
    E_time_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        # time='2012-02-13T21:07:02.77Z'
        print(time)
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        dist_km = dat.iloc[i]['dist(km)']
        distance_in_degree = dist_km / 111
        B_time = a_abstime(ot, depth, distance_in_degree)
        # calculate Etime
        date, t = time.split('T')
        event_folder = ''.join(date.split('-')) + ''.join(t.split(':'))[:-1]
        st = obspy.read(os.path.join(data_path, event_folder, '*Z'))
        trace = st[0]
        sta, lta, thr_on, thr_off = sta_lta_parameter
        if fig_outpath is not None:
            fig_outfile = os.path.join(fig_outpath, event_folder + '.pdf')
        else:
            fig_outfile = None
        E_time = make_stalta_as_Etime_abs(
            trace, ot, sta, lta, thr_on, thr_off, fig_outfile)

        B_time_list.append(str(B_time))
        E_time_list.append(str(E_time))
    dat['B_time'] = B_time_list
    dat['E_time'] = E_time_list

    dat.to_csv(catalog_outfile, index=False)


def Bvel_Efixedwindow(catalogfile, timewindow, catalog_outfile):
    dat = pd.read_csv(catalogfile)
    # calculate B_time and E_time
    E_time_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        dist_km = dat.iloc[i]['dist(km)']
        distance_in_degree = dist_km / 111
        B_time = UTCDateTime(dat.iloc[i]['B_time'])
        E_time = B_time + timewindow

        E_time_list.append(E_time)
    dat['E_time'] = E_time_list

    dat.to_csv(catalog_outfile, index=False)
    return None


# Main functions for background catalogs.
def preprocess_Xiaojiang_absBEtime(dat_path, dat_file, timeWindow_len):
    dat = pd.read_csv(os.path.join(dat_path, dat_file), encoding="gb2312")
    # calculate B_time and E_time
    B_time_list = []
    E_time_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        dist_km = dat.iloc[i]['dist(km)']
        distance_in_degree = dist_km / 111
        B_time = a_abstime(ot, depth, distance_in_degree)

        E_time = B_time + timeWindow_len

        B_time_list.append(str(B_time))
        E_time_list.append(str(E_time))
    dat['B_time'] = B_time_list
    dat['E_time'] = E_time_list

    dat.to_csv(os.path.join(dat_path, dat_file), index=False)


def catalog_during_days(teleseismic_catalog, out_file, dayWindow=30):
    out_tele = []

    catalog = pd.read_csv(teleseismic_catalog)

    tele_ot = catalog['time'].values
    A_time = catalog['A_time'].values
    Tb_B = catalog['Tb_B_time'].values
    Tb_E = catalog['Tb_E_time'].values
    Te_B = catalog['Te_B_time'].values
    Te_E = catalog['Te_E_time'].values

    lat = catalog['latitude'].values
    lon = catalog['longitude'].values
    depth = catalog['depth'].values
    dist = catalog['dist(km)'].values
    f_min = catalog['f_min'].values
    f_max = catalog['f_max'].values
    for i in range(len(tele_ot)):
        tele_datetime = UTCDateTime(tele_ot[i])
        print ('Background days for: '+str(tele_datetime))
        A_time_datetime = UTCDateTime(A_time[i])
        Tb_B_datetime = UTCDateTime(Tb_B[i])
        Tb_E_datetime = UTCDateTime(Tb_E[i])
        Te_B_datetime = UTCDateTime(Te_B[i])
        Te_E_datetime = UTCDateTime(Te_E[i])

        days = np.arange(-1 * (dayWindow/2), (dayWindow/2) + 1, 1)
        for day in days[days != 0]:
            tar_time = tele_datetime + timedelta(days=int(day))
            tar_A_time = A_time_datetime + timedelta(days=int(day))
            tar_Tb_B_time = Tb_B_datetime + timedelta(days=int(day))
            tar_Tb_E_time = Tb_E_datetime + timedelta(days=int(day))
            tar_Te_B_time = Te_B_datetime + timedelta(days=int(day))
            tar_Te_E_time = Te_E_datetime + timedelta(days=int(day))
            out_tele.append([str(tar_time), str(tar_A_time), str(tar_Tb_B_time), str(
                tar_Tb_E_time), str(tar_Te_B_time), str(tar_Te_E_time), lat[i], lon[i], depth[i], dist[i], f_min[i], f_max[i]])
    out_dataframe = pd.DataFrame(
        data=out_tele,
        columns=[
            'time',
            'A_time',
            'Tb_B_time',
            'Tb_E_time',
            'Te_B_time',
            'Te_E_time',
            'latitude',
            'longitude',
            'depth',
            'dist(km)',
            'f_min',
            'f_max'])
    out_dataframe.to_csv(out_file, index=False)
    return None


if __name__ == '__main__':
    # #referenceLocation = [26.5, 103]
    # dat_path = '/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Xiaojiang'
    # dat_file = 'EQ_info.csv'
    # timeWindow_len=8000
    # #out_file = 'EQ_info_part2.csv'
    #
    # preprocess_Geysers_absBEtime(dat_path, dat_file)
    # #preprocess_Xiaojiang_absBEtime(dat_path, dat_file, timeWindow_len)
    catalogfile = '/Users/yunnaidan/Project/Research/Dynamic_triggering/Xiaojiang/Data/Geysers/vel5_timewindow/EQ_info_GeysersGDX.csv'
    #preprocess_Geysers_changeTimewindow(catalogfile)
    print('Finish')
