#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2020/04/12
@file: utils.py
"""
import re
import os
import numpy as np
import pandas as pd
import obspy as ob
from scipy.signal import welch, spectrogram
from scipy.integrate import simps
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from datetime import timedelta
from math import radians, cos, sin, asin, sqrt, ceil
import warnings

from seispy import waveform as wf


def gen_target_days(days, origin_times):
    days_before, days_after = days
    target_days_all = []
    for ot in origin_times:
        ot_date = UTCDateTime(ot.year, ot.month, ot.day)
        target_days = [
            str(ot_date + timedelta(days=int(day)))
            for day in np.arange(-1 * days_before, days_after + 1, 1)
        ]
        target_days_all += target_days

    target_days_unique = sorted(list(set(target_days_all)))
    target_days_date = np.array([UTCDateTime(day) for day in target_days_unique])
    return target_days_date


def haversine(lon1, lat1, lon2, lat2):  # degree
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Translate degree to radian.
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula.
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of the Earth.
    return c * r, c * r / 111


def arrival(ot, depth, distance_in_degree):
    model = TauPyModel(model="iasp91")
    arrivals = model.get_travel_times(source_depth_in_km=depth,
                                      distance_in_degree=distance_in_degree)
    arrival_time = ot + arrivals[0].time
    return arrival_time


def phase(ot, vel, distance_in_km):
    travel_time = distance_in_km / vel
    phase_time = ot + travel_time
    return phase_time


def gen_time_windows(
        catalog_file=None,
        sta_file=None,
        tb=18000,
        te_b_vel=5.0,
        te_e_vel=2.0,
        f_min=15,
        f_max=40,
        f_bin=5,
        out_file=None
):
    if catalog_file is None:
        raise ValueError('Please input the catalog file of remote earthquakes!')
    if sta_file is None:
        raise ValueError('Please input the station name!')
    if out_file is None:
        raise ValueError('Please input the output file!')

    catalog = pd.read_csv(catalog_file)
    sta_df = pd.read_csv(sta_file)

    out_df = pd.DataFrame(columns=[
        'time', 'sta_name', 'dist_km', 'Tb_Begin', 'Tb_End', 'Te_Begin', 'Te_End',
    ])
    for sta_i, sta in sta_df.iterrows():
        reference_lon = sta['lon']
        reference_lat = sta['lat']
        sta_name = '.'.join([sta['net'], sta['sta'], str(sta['loc']).zfill(2), sta['cha']])
        print(sta_name)
        dist = [haversine(reference_lon,
                          reference_lat,
                          catalog.iloc[i]['longitude'],
                          catalog.iloc[i]['latitude'])
                for i in range(len(catalog))]
        dist_km = np.array(dist)[:, 0]
        dist_degree = np.array(dist)[:, 1]

        dist_km_list = []
        time_list = []
        tb_b_list = []
        tb_e_list = []
        te_b_list = []
        te_e_list = []
        for row_index, row in catalog.iterrows():
            time = row['time']
            print(time)

            if dist_km[row_index] > row['range']:
                ot = UTCDateTime(time)
                depth = row['depth']

                a_time = arrival(ot, depth, dist_degree[row_index])
                tb_b = a_time - tb
                tb_e = a_time
                te_b = phase(ot, te_b_vel, dist_km[row_index])
                te_e = phase(ot, te_e_vel, dist_km[row_index])

                dist_km_list.append(str(dist_km[row_index]))
                time_list.append(str(ot))
                tb_b_list.append(str(tb_b))
                tb_e_list.append(str(tb_e))
                te_b_list.append(str(te_b))
                te_e_list.append(str(te_e))
            else:
                print('The distance is too smaller: %f' % dist_km[row_index])

        if len(time_list) > 0:
            tmp_df = pd.DataFrame(columns=['time'],
                                  data=time_list)
            tmp_df['Tb_Begin'] = tb_b_list
            tmp_df['Tb_End'] = tb_e_list
            tmp_df['Te_Begin'] = te_b_list
            tmp_df['Te_End'] = te_e_list
            tmp_df['dist_km'] = dist_km_list
            tmp_df['sta_name'] = [sta_name] * len(tmp_df)

            out_df = out_df.append(tmp_df)

    out_df['f_min'] = np.full(len(out_df), f_min)
    out_df['f_max'] = np.full(len(out_df), f_max)
    out_df['f_bin'] = np.full(len(out_df), f_bin)
    out_df.to_csv(out_file, index=False)

    return None


def load_gf(sac_file, gf_info_file):
    df = pd.read_csv(gf_info_file)
    pz_file = df.loc[(df['sacfile'].values == sac_file), 'PZ_file'].values
    if len(pz_file) == 0 or str(pz_file[0]) == 'nan':
        warnings.warn('No sacPZ file for %s!' % sac_file)
        record_type = None
        sensitivity = None
        normalizing = None
        zeros = None
        poles = None
    else:
        pz_file = pz_file[0]
        with open(pz_file, 'r') as f:
            text = f.readlines()

        type_line = text[21]
        record_type = type_line.split(':')[1].split('(')[1].split(')')[0]

        sensitivity_line = text[21]
        sensitivity = float(sensitivity_line.split(':')[1].split('(')[0])
        normalizing_line = text[22]
        normalizing = float(normalizing_line.split(':')[1][:-1])

        zero_no = int(re.split('\s+', text[24])[1])
        zeros = []
        for i in range(zero_no):
            zero_info = text[25 + i]
            real, im = list(filter(None, re.split('\s+', zero_info)))
            zeros.append(complex(float(real), float(im)))

        # Delete zero points equalling to zero according to the data type.
        if record_type == 'M/S':
            zeros.remove(0.0)
        if record_type == 'M/S**2':
            zeros.remove(0.0)
            zeros.remove(0.0)

        pole_line_index = 24 + zero_no + 1
        pole_no = int(re.split('\s+', text[pole_line_index])[1])
        poles = []
        for i in range(pole_no):
            pole_info = text[pole_line_index + 1 + i]
            real, im = list(filter(None, re.split('\s+', pole_info)))
            poles.append(complex(float(real), float(im)))

    return record_type, sensitivity, normalizing, zeros, poles


def gf(sensitivity, normalizing, zeros, poles, f):
    s = complex(0, 2 * np.pi * f)
    gf1 = 1
    for zero in zeros:
        gf1 = gf1 * (s - zero)
    gf2 = 1
    for pole in poles:
        gf2 = gf2 * (s - pole)
    gf = sensitivity * normalizing * gf1 / gf2

    return abs(gf)


def catalog_during_days(teleseismic_catalog, out_file, day_window):
    out_tele = []

    catalog = pd.read_csv(teleseismic_catalog)

    tele_ot = catalog['time'].values
    sta_name = catalog['sta_name'].values
    tb_b = catalog['Tb_Begin'].values
    tb_e = catalog['Tb_End'].values
    te_b = catalog['Te_Begin'].values
    te_e = catalog['Te_End'].values

    f_min = catalog['f_min'].values
    f_max = catalog['f_max'].values
    f_bin = catalog['f_bin'].values
    for i in range(len(tele_ot)):
        tele_datetime = UTCDateTime(tele_ot[i])

        tb_b_datetime = UTCDateTime(tb_b[i])
        tb_e_datetime = UTCDateTime(tb_e[i])
        te_b_datetime = UTCDateTime(te_b[i])
        te_e_datetime = UTCDateTime(te_e[i])

        days = np.arange(-1 * day_window[0], day_window[1] + 1, 1)
        for day in days[days != 0]:
            tar_time = tele_datetime + timedelta(days=int(day))
            tar_tb_b_time = tb_b_datetime + timedelta(days=int(day))
            tar_tb_e_time = tb_e_datetime + timedelta(days=int(day))
            tar_te_b_time = te_b_datetime + timedelta(days=int(day))
            tar_te_e_time = te_e_datetime + timedelta(days=int(day))
            out_tele.append([str(tar_time),
                             sta_name[i],
                             f_min[i],
                             f_max[i],
                             f_bin[i],
                             str(tar_tb_b_time),
                             str(tar_tb_e_time),
                             str(tar_te_b_time),
                             str(tar_te_e_time)
                             ])
    out_dataframe = pd.DataFrame(
        data=out_tele,
        columns=[
            'time',
            'sta_name',
            'f_min',
            'f_max',
            'f_bin',
            'Tb_Begin',
            'Tb_End',
            'Te_Begin',
            'Te_End'
        ])
    out_dataframe.to_csv(out_file, index=False)
    return None


def abs_time(day_date, time_segment, i):
    year = int(str(day_date)[:4])
    month = int(str(day_date)[4:6])
    day = int(str(day_date)[6:8])
    hour = int((i * time_segment) / 3600)
    minute = int(((i * time_segment) % 3600) / 60)
    second = int(((i * time_segment) % 3600) % 60)
    abs_time_value = UTCDateTime(year, month, day, hour, minute, second)
    return abs_time_value


def load_pi(full_sta, tb_b, te_e, pi_database_path):
    b_day = ''.join([str(tb_b.year),
                     str(tb_b.month).zfill(2),
                     str(tb_b.day).zfill(2)])
    e_day = ''.join([str(te_e.year),
                     str(te_e.month).zfill(2),
                     str(te_e.day).zfill(2)])
    file_exist = False
    if b_day == e_day:
        tar_file = os.path.join(
            pi_database_path,
            'PI_' +
            full_sta +
            '_' +
            str(b_day) +
            '.csv')
        if os.path.exists(tar_file):
            file_exist = True
            pi_data_frame = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    'PI_' +
                    full_sta +
                    '_' +
                    str(b_day) +
                    '.csv'))
        else:
            pi_data_frame = None
    else:
        tar_file1 = os.path.join(
            pi_database_path, full_sta + '_' + str(b_day) + '.csv')
        tar_file2 = os.path.join(
            pi_database_path, full_sta + '_' + str(e_day) + '.csv')
        if os.path.exists(tar_file1) and os.path.exists(tar_file2):
            file_exist = True
            pi_data_frame1 = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    full_sta +
                    '_' +
                    str(b_day) +
                    '.csv'))
            pi_data_frame2 = pd.read_csv(
                os.path.join(
                    pi_database_path,
                    full_sta +
                    '_' +
                    str(e_day) +
                    '.csv'))
            pi_data_frame = pi_data_frame1.append(
                pi_data_frame2, ignore_index=True)
        else:
            pi_data_frame = None
    return file_exist, pi_data_frame


def select_event_waveforms(full_sta, tb_b, te_e, waveform_path, return_sac_file=False):
    b_day = ''.join([str(tb_b.year),
                     str(tb_b.month).zfill(2),
                     str(tb_b.day).zfill(2)])
    e_day = ''.join([str(te_e.year),
                     str(te_e.month).zfill(2),
                     str(te_e.day).zfill(2)])

    e_day_waveform = os.path.join(waveform_path,
                                  str(te_e.year),
                                  e_day,
                                  e_day + '.' + full_sta)

    if not os.path.exists(e_day_waveform):
        warnings.warn('No event waveform for %s!' % full_sta)
        # event_tr = None
        event_st = None
        sac_file = None
    else:
        if b_day == e_day:
            event_st = ob.read(e_day_waveform)
            sac_file = [e_day_waveform]
            # event_tr = event_st[0]
        else:
            b_day_waveform = os.path.join(waveform_path,
                                          str(tb_b.year),
                                          b_day,
                                          b_day + '.' + full_sta)
            if not os.path.exists(b_day_waveform):
                warnings.warn('No waveform of the day before the event for %s!' % full_sta)
                # event_tr = None
                event_st = None
                sac_file = None
            else:
                event_st = ob.read(b_day_waveform)
                event_st += ob.read(e_day_waveform)
                # event_st.merge(method=0, fill_value=0.0)
                # event_tr = event_st[0]
                sac_file = [b_day_waveform, e_day_waveform]
    if return_sac_file:
        output = [event_st, sac_file]
    else:
        output = event_st
    return output


def _gen_f_win_list(
        f_min, f_step, f_max
):
    f_win_no = int((f_max-f_min)/f_step)
    f_win_list = np.zeros([f_win_no, 2])
    for f_win_index in range(f_win_no):
        f_win_list[f_win_index] = [f_win_index*f_step + f_min,
                                   (f_win_index+1)*f_step + f_min]

    return f_win_list


def _f_column_name(
        f_min, f_step, f_max
):
    f_win_list = _gen_f_win_list(f_min, f_step, f_max)
    column_name = []
    for f_win in f_win_list:
        name = str(int(f_win[0])) + '-' + str(int(f_win[1]))
        column_name.append(name)
    return column_name


'''
def _tar_segments(st, tb_b, tb_e, te_b, te_e, time_segment):
    st_date = st
    et_date = te_e.date
    total_time = (et_date - st_date).seconds + 86400
    segment_count = int(total_time / time_segment)
    segment_times = np.array([abs_time(st.strftime('%Y%m%d'), time_segment, i)
                              for i in range(segment_count)])
    tar_segments_b0 = np.where((segment_times >= tb_b) &
                               (segment_times <= tb_e))[0]
    tar_segments_b = segment_times[tar_segments_b0 - 1]
    tar_segments_e = segment_times[(segment_times >= te_b) &
                                   (segment_times <= te_e)]

    return tar_segments_b, tar_segments_e


def _segments_psd(tr, segments, time_segment, nfft, gf_paras):
    sensitivity, normalizing, zeros, poles = gf_paras
    pxx_array = []
    for seg in segments:
        tar_tr = tr.slice(seg, seg + time_segment)
        data = tar_tr.data
        fs = tar_tr.stats.sampling_rate
        pxx_all, f_all = psd(data, fs, nfft)

        gf_list = np.array([gf(sensitivity, normalizing, zeros, poles, f) for f in f_all])
        ground_pxx = np.array([(pxx_all[i] / (gf_list[i] ** 2)) * (10 ** 18) for i in range(len(gf_list))])

        pxx_array.append(ground_pxx)

    pxx_array = np.array(pxx_array)
    pxx_mean = np.mean(pxx_array, axis=0)

    return f_all, pxx_mean


def window_psds(
        full_sta=None,
        tb_b=None,
        tb_e=None,
        te_b=None,
        te_e=None,
        waveform_path=None,
        time_segment=None,
        nfft=None,
        gf_paras=None,
):
    tr = wf.select_event_waveform(full_sta, tb_b, te_e, waveform_path)
    tr.detrend('linear')
    tr.detrend('constant')

    st_date = tr.stats.starttime.date

    tar_segments_b, tar_segments_e = \
        _tar_segments(st_date, tb_b, tb_e, te_b, te_e, time_segment)
    f_b, ground_pxx_mean_b = _segments_psd(tr, tar_segments_b, time_segment, nfft, gf_paras)
    f_e, ground_pxx_mean_e = _segments_psd(tr, tar_segments_e, time_segment, nfft, gf_paras)

    return f_b, ground_pxx_mean_b, f_e, ground_pxx_mean_e
'''


def get_act_windows(
        event=None,
        full_sta=None,
        pi_database_folder=None,
        tb_b=None,
        tb_e=None,
        te_b=None,
        te_e=None,
        sperseg=None,
):
    file_exist, pi_data_frame = load_pi(
        full_sta, tb_b, te_e, pi_database_folder)
    if file_exist:
        pi_data_frame['time'] = [UTCDateTime(t) for t in pi_data_frame['time']]
        tar_df_b = pi_data_frame[(pi_data_frame['time'] >= tb_b) & (pi_data_frame['time'] <= tb_e)]
        tar_df_e = pi_data_frame[(pi_data_frame['time'] >= te_b) & (
                pi_data_frame['time'] <= te_e)]
        act_tb_b = tar_df_b.iloc[0]['time'] - sperseg / 2
        act_tb_e = tar_df_b.iloc[-1]['time'] + sperseg / 2
        act_te_b = tar_df_e.iloc[0]['time'] - sperseg / 2
        act_te_e = tar_df_e.iloc[-1]['time'] + sperseg / 2

    else:
        warnings.warn('No PI data of %s for %s.' % (full_sta, event))
        act_tb_b = None
        act_tb_e = None
        act_te_b = None
        act_te_e = None
    return act_tb_b, act_tb_e, act_te_b, act_te_e


def high_frequency_energy(
        f_b=None,
        ground_pxx_mean_b=None,
        f_e=None,
        ground_pxx_mean_e=None,
        f_min=None,
        f_max=None,
        f_step=None,
):
    f_win_list = _gen_f_win_list(f_min, f_step, f_max)

    pi_b = []
    pi_e = []
    for f_win in f_win_list:
        f_min, f_max = f_win

        tar_index_b = np.where((f_b >= f_min) & (f_b <= f_max))[0]
        tar_index_e = np.where((f_e >= f_min) & (f_e <= f_max))[0]

        pi_b.append(simps(ground_pxx_mean_b[tar_index_b], f_b[tar_index_b]))
        pi_e.append(simps(ground_pxx_mean_e[tar_index_e], f_e[tar_index_e]))

    high_pi_b = np.sum(pi_b)
    high_pi_e = np.sum(pi_e)

    return high_pi_b, high_pi_e


