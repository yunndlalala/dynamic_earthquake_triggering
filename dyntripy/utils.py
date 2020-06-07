import re
import numpy as np
import pandas as pd
from scipy.signal import welch
from obspy import UTCDateTime
from datetime import timedelta
from math import radians, cos, sin, asin, sqrt, ceil


def gen_target_days(days, origin_times):
    days_before, days_after = days
    target_days_all = []
    for ot in origin_times:
        ot_date = UTCDateTime(ot.year, ot.month, ot.day)
        target_days = [str(ot_date + timedelta(days=int(day)))
                           for day in np.arange(-1 * days_before, days_after + 1, 1)]
        target_days_all += target_days

    target_days_unique = sorted(list(set(target_days_all)))
    target_days_date = np.array([UTCDateTime(day) for day in target_days_unique])
    return target_days_date


def psd(data, fs):
    nfft = 512
    seg = ceil(len(data) / (nfft / 2))
    nperseg = int(len(data) / seg) * 2
    f, pxx_all = welch(data,
                       fs,
                       window='hanning',
                       nperseg=nperseg,
                       noverlap=nperseg / 2,
                       nfft=nfft,
                       detrend=None,
                       return_onesided=True,
                       scaling='density',
                       axis=-1)

    return pxx_all, f


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


def load_gf(sac_file, gf_info_file):
    df = pd.read_csv(gf_info_file)
    PZ_file = df.loc[(df['sacfile'].values == sac_file), 'PZ_file'].values[0]

    with open(PZ_file, 'r') as f:
        text = f.readlines()

    type_line = text[21]
    type = type_line.split(':')[1].split('(')[1].split(')')[0]

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
    if type == 'M/S':
        zeros.remove(0.0)
    if type == 'M/S**2':
        zeros.remove(0.0)
        zeros.remove(0.0)

    pole_line_index = 24 + zero_no + 1
    pole_no = int(re.split('\s+', text[pole_line_index])[1])
    poles = []
    for i in range(pole_no):
        pole_info = text[pole_line_index + 1 + i]
        real, im = list(filter(None, re.split('\s+', pole_info)))
        poles.append(complex(float(real), float(im)))

    return type, sensitivity, normalizing, zeros, poles


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


def catalog_during_days(teleseismic_catalog, out_file, day_window=[60, 60]):
    out_tele = []

    catalog = pd.read_csv(teleseismic_catalog)

    tele_ot = catalog['time'].values
    tb_b = catalog['Tb_Begin'].values
    tb_e = catalog['Tb_End'].values
    te_b = catalog['Te_Begin'].values
    te_e = catalog['Te_End'].values

    f_min = catalog['f_min'].values
    f_max = catalog['f_max'].values
    for i in range(len(tele_ot)):
        tele_datetime = UTCDateTime(tele_ot[i])
        # print('Background days for: ' + str(tele_datetime), end='\r')
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
                             f_min[i],
                             f_max[i],
                             str(tar_tb_b_time),
                             str(tar_tb_e_time),
                             str(tar_te_b_time),
                             str(tar_te_e_time)
                             ])
    out_dataframe = pd.DataFrame(
        data=out_tele,
        columns=[
            'time',
            'f_min',
            'f_max',
            'Tb_Begin',
            'Tb_End',
            'Te_Begin',
            'Te_End'
        ])
    out_dataframe.to_csv(out_file, index=False)
    return None
