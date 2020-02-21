import re
import numpy as np
import pandas as pd
import obspy
from scipy.signal import welch
from math import radians, cos, sin, asin, sqrt, ceil


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


def peak_dynamic_stress(tele_event_tr, shear_wave_velocity, te_b, te_e, shear_modulus=3e10):
    tr = tele_event_tr.copy()

    fs = tr.stats.sampling_rate
    abs_data_begin = tr.stats.starttime
    te_b_sec = float('%.2f' % (te_b - abs_data_begin))
    te_e_sec = float('%.2f' % (te_e - abs_data_begin))

    data = tr[round(te_b_sec*fs):round(te_e_sec*fs)]
    peak_ground_record = np.max(data)

    dynamic_stress = 1e-12 * (shear_modulus * peak_ground_record) / shear_wave_velocity

    return dynamic_stress


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
        _, real, im = re.split('\s+', zero_info[:-1])
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
        _, real, im = re.split('\s+', pole_info[:-1])
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


if __name__ == "__main__":
    print(haversine(103, 26.5, -76.36, 1.93))
    print('finish')
