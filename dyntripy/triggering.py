#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: run_utils.py
"""
import os
import re
import numpy as np
import pandas as pd
from obspy import UTCDateTime

from dyntripy.hifi import PIdatabase, HighPIRatio, AssociateBackgroundPIR, ConfidenceLevel
from dyntripy.utils import gen_target_days, catalog_during_days


class Triggering(object):

    def __init__(self, hypers):
        self.hypers = hypers
        self.sta_list = self.load_sta()
        self.f_win_list = self.gen_f_win_list()

    def load_sta(self):
        with open(self.hypers['data_source']['station_file'], 'r') as f:
            text = f.readlines()
        sta_list = []
        for line in text[1:]:
            line = re.split('\s+', line)[0]
            net, sta, loc, chn = line.split(',')
            full_sta = '.'.join([net, sta, loc, chn])
            if full_sta not in sta_list:
                sta_list.append(full_sta)

        return sta_list

    def gen_f_win_list(self):
        f_min, f_step, f_max = self.hypers['net_database']['frequency_segment']

        f_win_no = int((f_max-f_min)/f_step)
        f_win_list = np.zeros([f_win_no, 2])
        for f_win_index in range(f_win_no):
            f_win_list[f_win_index] = [f_win_index*f_step + f_min,
                                       (f_win_index+1)*f_step + f_min]

        return f_win_list

    def net_database(self, p=1):
        if not os.path.exists(self.hypers['net_database']['database_path']):
            os.makedirs(self.hypers['net_database']['database_path'])

        catalog_df = pd.read_csv(self.hypers['data_source']['remote_earthquake_catalog'])
        origin_times = np.array([UTCDateTime(t) for t in catalog_df['time'].values])
        target_days = gen_target_days(self.hypers['net_database']['days'], origin_times)

        for full_sta in self.sta_list:
            print (full_sta)
            out_folder = os.path.join(self.hypers['net_database']['database_path'], full_sta)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            PIdatabase.run_pi_parallel(
                full_sta,
                target_days,
                self.hypers['data_source']['waveform_path'],
                self.hypers['net_database']['time_segment'],
                self.f_win_list,
                self.hypers['data_source']['response_file'],
                out_folder,
                p)
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
                print ('Build background catalog ...\n')
                self.bg_catalog()
            out_path = self.hypers['net_ratio']['RB_path']
            out_file_prefix = 'RB'

        if not os.path.exists(out_path):
            os.makedirs(out_path)

        for full_sta in self.sta_list:
            print(full_sta)
            pi_database_folder = os.path.join(self.hypers['net_database']['database_path'],
                                              full_sta)
            out_file = os.path.join(out_path, out_file_prefix + '_' + full_sta + '.csv')
            HighPIRatio.run_pir_parallel(
                self.hypers['net_database']['frequency_segment'][1],
                catalog,
                pi_database_folder,
                full_sta,
                out_file,
                p)
        out_df = pd.read_csv(out_file)
        out_df.sort_values(by='time', inplace=True)
        out_df.to_csv(out_file, index=False)

        return None

    def match_ratio(self, full_sta, p=1):
        tele_pir_file = os.path.join(
            self.hypers['net_ratio']['RE_path'], 'RE_' + full_sta + '.csv')
        background_pir_file = os.path.join(
            self.hypers['net_ratio']['RB_path'], 'RB_' + full_sta + '.csv')
        out_file = os.path.join(
            self.hypers['net_cl']['matched_ratio_path'], 'matched_' + full_sta + '.csv')
        AssociateBackgroundPIR.run_associate_parallel(
            self.hypers['data_source']['remote_earthquake_catalog'],
            tele_pir_file,
            background_pir_file,
            'ratio',
            self.hypers['net_ratio']['background_days'],
            out_file,
            p)

        return None

    def net_cl(self, p=1):
        if not os.path.exists(self.hypers['net_cl']['cl_path']):
            os.makedirs(self.hypers['net_cl']['cl_path'])
        if not os.path.exists(self.hypers['net_cl']['matched_ratio_path']):
            os.makedirs(self.hypers['net_cl']['matched_ratio_path'])
        if self.hypers['net_cl']['figure_path'] != 'None':
            if not os.path.exists(self.hypers['net_cl']['figure_path']):
                os.makedirs(self.hypers['net_cl']['figure_path'])

        for full_sta in self.sta_list:
            print (full_sta)
            background_pir_associated_file = os.path.join(
                self.hypers['net_cl']['matched_ratio_path'],
                'matched_' + full_sta + '.csv')
            if not os.path.exists(background_pir_associated_file):
                print ('Match ratios ...\n')
                self.match_ratio(full_sta, p=p)
                matched_df = pd.read_csv(background_pir_associated_file)
                matched_df.sort_values(by='time', inplace=True)
                matched_df.to_csv(background_pir_associated_file, index=False)

            out_file = os.path.join(
                self.hypers['net_cl']['cl_path'], 'CL_' + full_sta + '.csv')

            if self.hypers['net_cl']['figure_path'] != 'None':
                figure_out_folder = os.path.join(self.hypers['net_cl']['figure_path'], full_sta)
                if not os.path.exists(figure_out_folder):
                    os.makedirs(figure_out_folder)
            else:
                figure_out_folder = None

            ConfidenceLevel.run_cl_parallel(
                background_pir_associated_file,
                self.hypers['net_cl']['threshold'],
                out_file,
                figure_out_folder,
                p)
            cl_df = pd.read_csv(out_file)
            cl_df.sort_values(by='time', inplace=True)
            cl_df.to_csv(out_file, index=False)

        return None


