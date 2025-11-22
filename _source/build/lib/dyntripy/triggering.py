#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: run_utils.py
"""
import os
import numpy as np
import pandas as pd
from obspy import UTCDateTime

from dyntripy.hifi import PIdatabase, HighPIRatio, AssociateBackgroundPIR, ConfidenceLevel
from dyntripy.utils import gen_target_days, catalog_during_days, _gen_f_win_list


class Triggering(object):

    def __init__(self, hypers):
        self.hypers = hypers
        self.sta_list = self.load_sta()

    def load_sta(self):
        sta_df = pd.read_csv(self.hypers['data_source']['station_file'], dtype=str)

        sta_list = np.array([
            '.'.join([row['net'], row['sta'], row['loc'], row['cha']])
            for row_i, row in sta_df.iterrows()
        ])

        return sta_list

    def net_database(self, p=1):
        if not os.path.exists(self.hypers['net_database']['database_path']):
            os.makedirs(self.hypers['net_database']['database_path'])

        catalog_df = pd.read_csv(self.hypers['data_source']['remote_earthquake_catalog'])

        f_min, f_step, f_max = self.hypers['net_database']['frequency_segment']
        f_win_list = _gen_f_win_list(f_min, f_step, f_max)

        if len(str(f_min).split('.')) > 1 \
                or len(str(f_max).split('.')) > 1 \
                or len(str(f_step).split('.')) > 1:
            raise ValueError('The min, max and step values of frequency must be integer.')

        for full_sta in self.sta_list:
            print(full_sta)

            tar_catalog_df = catalog_df[catalog_df['sta_name'] == full_sta]
            origin_times = np.array([UTCDateTime(t) for t in tar_catalog_df['time'].values])
            target_days = gen_target_days(self.hypers['net_database']['days'], origin_times)
            out_folder = os.path.join(self.hypers['net_database']['database_path'], full_sta)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            PIdatabase.run_pi_parallel(
                full_sta,
                target_days,
                self.hypers['data_source']['waveform_path'],
                self.hypers['net_database']['nperseg'],
                self.hypers['net_database']['noverlap'],
                self.hypers['net_database']['nfft'],
                f_win_list,
                self.hypers['data_source']['response_file'],
                out_folder,
                p
            )

        return None

    def bg_catalog(self):
        catalog_during_days(self.hypers['data_source']['remote_earthquake_catalog'],
                            self.hypers['net_ratio']['background_catalog'],
                            self.hypers['net_ratio']['background_days'])
        return None

    def net_ratio(self, bg=False, p=1):
        if not bg:
            catalog = self.hypers['data_source']['remote_earthquake_catalog']
            out_path = self.hypers['net_ratio']['RE_path']
            out_file_prefix = 'RE'

        else:
            catalog = self.hypers['net_ratio']['background_catalog']
            if not os.path.exists(catalog):
                print('Build background catalog ...\n')
                self.bg_catalog()
            out_path = self.hypers['net_ratio']['RB_path']
            out_file_prefix = 'RB'

        if not os.path.exists(out_path):
            os.makedirs(out_path)

        HighPIRatio.run_pir_parallel(
            self.hypers['net_database']['frequency_segment'][1],
            catalog,
            self.hypers['net_database']['database_path'],
            out_path,
            out_file_prefix,
            p
        )

        return None

    def net_match(self, p=1):
        if not os.path.exists(self.hypers['net_cl']['matched_ratio_path']):
            os.makedirs(self.hypers['net_cl']['matched_ratio_path'])

        AssociateBackgroundPIR.run_associate_parallel(
            self.hypers['data_source']['remote_earthquake_catalog'],
            'ratio',
            self.hypers['net_ratio']['background_days'],
            self.hypers['net_ratio']['RE_path'],
            self.hypers['net_ratio']['RB_path'],
            self.hypers['net_cl']['matched_ratio_path'],
            p,
        )

        return None

    def net_cl(self, p=1):
        if not os.path.exists(self.hypers['net_cl']['cl_path']):
            os.makedirs(self.hypers['net_cl']['cl_path'])
        if self.hypers['net_cl']['figure_path'] != 'None':
            if not os.path.exists(self.hypers['net_cl']['figure_path']):
                os.makedirs(self.hypers['net_cl']['figure_path'])

        ConfidenceLevel.run_cl_parallel(
            self.hypers['data_source']['remote_earthquake_catalog'],
            self.hypers['net_cl']['matched_ratio_path'],
            self.hypers['net_cl']['cl_path'],
            self.hypers['net_cl']['figure_path'],
            self.hypers['net_cl']['threshold'],
            p
        )

        return None


