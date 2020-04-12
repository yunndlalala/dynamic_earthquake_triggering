import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from obspy.taup import TauPyModel

from dyntrigger.utils.basic_utils import haversine



# TOOL
def str2datetime(str_list):
    datetime_lst = [datetime.strptime(
        i[:10] + ' ' + i[11:-1], '%Y-%m-%d %H:%M:%S.%f') for i in str_list]
    return datetime_lst


def abs_arrival_time(ot, depth, distance_in_degree):
    model = TauPyModel(model="iasp91")
    arrivals = model.get_travel_times(source_depth_in_km=depth,
                                      distance_in_degree=distance_in_degree)
    arrival_time = ot + arrivals[0].time
    return arrival_time


def abs_phase_time(ot, vel, distance_in_km):
    travel_time = distance_in_km / vel
    phase_time = ot + travel_time
    return phase_time


# Main functions for time windows
def begin_a_end_fixed(catalog_file, sta_lat, sta_lon, tb, te):
    dat = pd.read_csv(catalog_file)

    a_time_list = []
    tb_b_list = []
    tb_e_list = []
    te_b_list = []
    te_e_list = []
    dist_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        event_lat = dat.iloc[i]['latitude']
        event_lon = dat.iloc[i]['longitude']
        if sta_lon is None and sta_lat is None:
            dist_km = dat.iloc[i]['dist(km)']
            distance_in_degree = dist_km / 111
        else:
            dist_km, distance_in_degree = haversine(
                sta_lon, sta_lat, event_lon, event_lat)

        a_time = abs_arrival_time(ot, depth, distance_in_degree)
        tb_b = a_time - tb
        tb_e = a_time
        te_b = a_time
        te_e = te_b + te

        a_time_list.append(str(a_time))
        tb_b_list.append(str(tb_b))
        tb_e_list.append(str(tb_e))
        te_b_list.append(str(te_b))
        te_e_list.append(str(te_e))
        dist_list.append(dist_km)

    dat['A_time'] = a_time_list
    dat['Tb_B_time'] = tb_b_list
    dat['Tb_E_time'] = tb_e_list
    dat['Te_B_time'] = te_b_list
    dat['Te_E_time'] = te_e_list
    if sta_lon is not None and sta_lat is not None:
        dat['dist(km)'] = dist_list

    dat.to_csv(catalog_file, index=False)

    return None



def begin_v_end_fixed(catalog_file, sta_lat, sta_lon, tb, vel_begin, te):
    dat = pd.read_csv(catalog_file)

    a_time_list = []
    tb_b_list = []
    tb_e_list = []
    te_b_list = []
    te_e_list = []
    dist_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        event_lat = dat.iloc[i]['latitude']
        event_lon = dat.iloc[i]['longitude']
        if sta_lon is None and sta_lat is None:
            dist_km = dat.iloc[i]['dist(km)']
            distance_in_degree = dist_km / 111
        else:
            dist_km, distance_in_degree = haversine(
                sta_lon, sta_lat, event_lon, event_lat)

        a_time = abs_arrival_time(ot, depth, distance_in_degree)
        tb_b = a_time - tb
        tb_e = a_time
        te_b = abs_phase_time(ot, vel_begin, dist_km)
        te_e = te_b + te

        a_time_list.append(str(a_time))
        tb_b_list.append(str(tb_b))
        tb_e_list.append(str(tb_e))
        te_b_list.append(str(te_b))
        te_e_list.append(str(te_e))
        dist_list.append(dist_km)

    dat['A_time'] = a_time_list
    dat['Tb_B_time'] = tb_b_list
    dat['Tb_E_time'] = tb_e_list
    dat['Te_B_time'] = te_b_list
    dat['Te_E_time'] = te_e_list
    if sta_lon is not None and sta_lat is not None:
        dat['dist(km)'] = dist_list

    dat.to_csv(catalog_file, index=False)

    return None


def begin_v_end_v(catalog_file, sta_lat, sta_lon, tb, vel_begin, vel_end):
    dat = pd.read_csv(catalog_file)

    a_time_list = []
    tb_b_list = []
    tb_e_list = []
    te_b_list = []
    te_e_list = []
    dist_list = []
    for i in range(len(dat)):
        time = dat.iloc[i]['time']
        ot = UTCDateTime(time)
        depth = dat.iloc[i]['depth']
        event_lat = dat.iloc[i]['latitude']
        event_lon = dat.iloc[i]['longitude']
        if sta_lon is None and sta_lat is None:
            dist_km = dat.iloc[i]['dist(km)']
            distance_in_degree = dist_km / 111
        else:
            dist_km, distance_in_degree = spherical_dist(
                sta_lon, sta_lat, event_lon, event_lat)

        a_time = abs_arrival_time(ot, depth, distance_in_degree)
        tb_b = a_time - tb
        tb_e = a_time
        te_b = abs_phase_time(ot, vel_begin, dist_km)
        te_e = abs_phase_time(ot, vel_end, dist_km)

        a_time_list.append(str(a_time))
        tb_b_list.append(str(tb_b))
        tb_e_list.append(str(tb_e))
        te_b_list.append(str(te_b))
        te_e_list.append(str(te_e))
        dist_list.append(dist_km)

    dat['A_time'] = a_time_list
    dat['Tb_B_time'] = tb_b_list
    dat['Tb_E_time'] = tb_e_list
    dat['Te_B_time'] = te_b_list
    dat['Te_E_time'] = te_e_list
    if sta_lon is not None and sta_lat is not None:
        dat['dist(km)'] = dist_list

    dat.to_csv(catalog_file, index=False)

    return None


# Main functions for background catalogs.
def catalog_during_days(teleseismic_catalog, out_file, day_window=30):
    out_tele = []

    catalog = pd.read_csv(teleseismic_catalog)

    tele_ot = catalog['time'].values
    a_time = catalog['A_time'].values
    tb_b = catalog['Tb_B_time'].values
    tb_e = catalog['Tb_E_time'].values
    te_b = catalog['Te_B_time'].values
    te_e = catalog['Te_E_time'].values

    lat = catalog['latitude'].values
    lon = catalog['longitude'].values
    depth = catalog['depth'].values
    dist = catalog['dist(km)'].values
    f_min = catalog['f_min'].values
    f_max = catalog['f_max'].values
    for i in range(len(tele_ot)):
        tele_datetime = UTCDateTime(tele_ot[i])
        print('Background days for: ' + str(tele_datetime))
        a_time_datetime = UTCDateTime(a_time[i])
        tb_b_datetime = UTCDateTime(tb_b[i])
        tb_e_datetime = UTCDateTime(tb_e[i])
        te_b_datetime = UTCDateTime(te_b[i])
        te_e_datetime = UTCDateTime(te_e[i])

        days = np.arange(-1 * (day_window / 2), (day_window / 2) + 1, 1)
        for day in days[days != 0]:
            tar_time = tele_datetime + timedelta(days=int(day))
            tar_a_time = a_time_datetime + timedelta(days=int(day))
            tar_tb_b_time = tb_b_datetime + timedelta(days=int(day))
            tar_tb_e_time = tb_e_datetime + timedelta(days=int(day))
            tar_te_b_time = te_b_datetime + timedelta(days=int(day))
            tar_te_e_time = te_e_datetime + timedelta(days=int(day))
            out_tele.append([str(tar_time),
                             str(tar_a_time),
                             str(tar_tb_b_time),
                             str(tar_tb_e_time),
                             str(tar_te_b_time),
                             str(tar_te_e_time),
                             lat[i],
                             lon[i],
                             depth[i],
                             dist[i],
                             f_min[i],
                             f_max[i]])
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
