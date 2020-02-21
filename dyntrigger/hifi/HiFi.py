#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: run_utils.py
"""
import os
import re
import time
import numpy as np
import pandas as pd
import obspy
import multiprocessing

from dyntrigger.hifi import PIdatabase, HighPIRatio, AssociateBackgroundPIR, ConfidenceLevel
from dyntrigger.utils.preprocess_tele_catalog import catalog_during_days
from dyntrigger.plot_utils.plot_hifi_result import plot_psd, plot_bg_pir
from dyntrigger.utils.basic_utils import peak_dynamic_stress


class HiFi(object):

    def __init__(self, hypers):
        self.hypers = hypers
        self.root_path = hypers['info']['root_path']

        self.sta_file = os.path.join(self.root_path, hypers['info']['sta_file'])
        self.sta_list = self.load_sta()

        self.chn_list = hypers['info']['chn_list']

        self.f_range = hypers['pi']['f_range']
        self.f_step = hypers['pi']['f_step']
        self.f_win_list = self.gen_f_win_list()

        self.time_segment = hypers['pi']['time_segment']
        self.gf_info = os.path.join(self.root_path, hypers['pi']['gf_info'])
        self.raw_data_path = os.path.join(self.root_path, hypers['pi']['raw_data_path'])
        self.pi_database_path = os.path.join(self.root_path, hypers['pi']['pi_database_path'])

        self.tele_catalog = os.path.join(self.root_path, hypers['tele_pir']['tele_catalog'])
        self.tele_pir_path = os.path.join(self.root_path, hypers['tele_pir']['tele_pir_path'])

        self.day_window = hypers['background_pir']['day_window']
        self.background_catalog = os.path.join(self.root_path, hypers['background_pir']['background_catalog'])
        self.background_pir_path = os.path.join(self.root_path, hypers['background_pir']['background_pir_path'])
        self.background_pir_associated_path = os.path.join(self.root_path,
                                                           hypers['background_pir']['background_pir_associated_path'])

        self.cl_path = os.path.join(self.root_path, hypers['confidence_level']['cl_path'])

        self.sta = hypers['plot']['sta']
        self.chn = hypers['plot']['chn']
        self.event_data_path_psd = os.path.join(self.root_path, hypers['plot']['event_data_path'])
        self.psd_figure_path = os.path.join(self.root_path, hypers['plot']['psd_figure_path'])
        self.background_pir_figure_path = os.path.join(self.root_path, hypers['plot']['background_pir_figure_path'])

        self.shear_wave_velocity = hypers['stress']['shear_wave_velocity']
        self.shear_modulus = hypers['stress']['shear_modulus']
        self.event_data_path_stress = os.path.join(self.root_path, hypers['stress']['event_data_path'])
        self.peak_dynamic_stress_file = os.path.join(self.root_path, hypers['stress']['peak_dynamic_stress_file'])

    def load_sta(self):
        with open(self.sta_file, 'r') as f:
            text = f.readlines()
        sta_list = []
        for line in text[1:]:
            net, sta, _, _, _ = line.split(',')
            full_sta = '.'.join([net, sta])
            if full_sta not in sta_list:
                sta_list.append(full_sta)

        return sta_list

    def gen_f_win_list(self):
        f_range = self.f_range
        f_step = self.f_step

        f_win_no = int((f_range[1]-f_range[0])/f_step)
        f_win_list = np.zeros([f_win_no, 2])
        for f_win_index in range(f_win_no):
            f_win_list[f_win_index] = [f_win_index*f_step + f_range[0],
                                       (f_win_index+1)*f_step + f_range[0]]

        return f_win_list

    def net_pi(self, cores=1):
        if os.path.exists(self.pi_database_path):
            os.removedirs(self.pi_database_path)
        os.makedirs(self.pi_database_path)

        for sta in self.sta_list:
            for chn in self.chn_list:
                print(sta, chn)
                out_folder = os.path.join(self.pi_database_path, sta + '.' + chn)
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)
                PIdatabase.run_pi_parallel(
                    sta,
                    chn,
                    self.raw_data_path,
                    self.time_segment,
                    self.f_win_list,
                    self.gf_info,
                    out_folder,
                    cores)
        return None

    def net_tele_pir(self, cores=1):
        if os.path.exists(self.tele_pir_path):
            os.removedirs(self.tele_pir_path)
        os.makedirs(self.tele_pir_path)

        if cores == 1:
            for sta in self.sta_list:
                for chn in self.chn_list:
                    print(sta, chn)
                    pi_database_folder = os.path.join(self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(self.tele_pir_path, sta + '.' + chn + '_PIR.csv')
                    HighPIRatio.run_pir(
                        self.f_step,
                        self.tele_catalog,
                        pi_database_folder,
                        sta,
                        chn,
                        out_file)
        else:
            pool = multiprocessing.Pool(processes=cores)
            tasks = []
            for sta in self.sta_list:
                for chn in self.chn_list:
                    pi_database_folder = os.path.join(self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(self.tele_pir_path, sta + '.' + chn + '_PIR.csv')
                    tasks.append(
                        (self.f_step,
                         self.tele_catalog,
                         pi_database_folder,
                         sta,
                         chn,
                         out_file))

            # chunksize is how many tasks will be processed by one processor
            rs = pool.starmap_async(HighPIRatio.run_pir, tasks, chunksize=1)
            # close() & join() is necessary
            pool.close()
            # simple progress bar
            while (True):
                remaining = rs._number_left
                print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
                      end='\r')  # '\r' means remove the last line
                if (rs.ready()):
                    break
                time.sleep(0.5)

            pool.join()

        return None

    def gen_background_catalog(self):
        catalog_during_days(self.tele_catalog, self.background_catalog, self.day_window)
        return None

    def net_background_pir(self,cores=1):
        if os.path.exists(self.background_pir_path):
            os.removedirs(self.background_pir_path)
        os.makedirs(self.background_pir_path)

        if cores == 1:
            for sta in self.sta_list:
                for chn in self.chn_list:
                    print(sta, chn)
                    pi_database_folder = os.path.join(
                        self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    HighPIRatio.run_pir(
                        self.f_step,
                        self.background_catalog,
                        pi_database_folder,
                        sta,
                        chn,
                        out_file)
        else:
            pool = multiprocessing.Pool(processes=cores)
            tasks = []
            for sta in self.sta_list:
                for chn in self.chn_list:
                    pi_database_folder = os.path.join(
                        self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    tasks.append(
                        (self.f_step,
                         self.background_catalog,
                         pi_database_folder,
                         sta,
                         chn,
                         out_file))

            # chunksize is how many tasks will be processed by one processor
            rs = pool.starmap_async(HighPIRatio.run_pir, tasks, chunksize=1)
            # close() & join() is necessary
            pool.close()
            # simple progress bar
            while (True):
                remaining = rs._number_left
                print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
                      end='\r')  # '\r' means remove the last line
                if (rs.ready()):
                    break
                time.sleep(0.5)

            pool.join()

        return None

    def net_background_pir_asso(self, cores=1):
        if os.path.exists(self.background_pir_associated_path):
            os.removedirs(self.background_pir_associated_path)
        os.makedirs(self.background_pir_associated_path)

        if cores == 1:
            for sta in self.sta_list:
                for chn in self.chn_list:
                    print(sta, chn)
                    background_pir_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    out_file = os.path.join(
                        self.background_pir_associated_path, sta + '.' + chn + '_background_PIR_associated.csv')
                    AssociateBackgroundPIR.run_associate(
                        self.tele_catalog, background_pir_file, 'PIR_log', self.day_window, out_file)
        else:
            pool = multiprocessing.Pool(processes=cores)
            tasks = []
            for sta in self.sta_list:
                for chn in self.chn_list:
                    background_pir_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    out_file = os.path.join(
                        self.background_pir_associated_path, sta + '.' + chn + '_background_PIR_associated.csv')
                    tasks.append(
                        (self.tele_catalog,
                         background_pir_file,
                         'PIR_log',
                         self.day_window,
                         out_file))

            # chunksize is how many tasks will be processed by one processor
            rs = pool.starmap_async(
                AssociateBackgroundPIR.run_associate,
                tasks,
                chunksize=1)
            # close() & join() is necessary
            pool.close()
            # simple progress bar
            while (True):
                remaining = rs._number_left
                print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
                      end='\r')  # '\r' means remove the last line
                if (rs.ready()):
                    break
                time.sleep(0.5)

            pool.join()

        return None

    def net_cl(self, cores=1):
        if os.path.exists(self.cl_path):
            os.removedirs(self.cl_path)
        os.makedirs(self.cl_path)

        if cores == 1:
            for sta in self.sta_list:
                for chn in self.chn_list:
                    print(sta, chn)
                    background_pir_associated_file = os.path.join(
                        self.background_pir_associated_path,
                        sta + '.' + chn + '_background_PIR_associated.csv')
                    tele_pir_file = os.path.join(
                        self.tele_pir_path, sta + '.' + chn + '_PIR.csv')
                    out_file = os.path.join(
                        self.cl_path, sta + '.' + chn + '_CL.csv')
                    ConfidenceLevel.run_cl(
                        background_pir_associated_file,
                        tele_pir_file,
                        out_file,
                        'PIR_log')
        else:
            pool = multiprocessing.Pool(processes=cores)
            tasks = []
            for sta in self.sta_list:
                for chn in self.chn_list:
                    background_pir_associated_file = os.path.join(
                        self.background_pir_associated_path,
                        sta + '.' + chn + '_background_PIR_associated.csv')
                    tele_pir_file = os.path.join(
                        self.tele_pir_path, sta + '.' + chn + '_PIR.csv')
                    out_file = os.path.join(
                        self.cl_path, sta + '.' + chn + '_CL.csv')
                    tasks.append((
                        background_pir_associated_file,
                        tele_pir_file,
                        out_file,
                        'PIR_log'))
            # chunksize is how many tasks will be processed by one processor
            rs = pool.starmap_async(
                ConfidenceLevel.run_cl,
                tasks,
                chunksize=1)
            # close() & join() is necessary
            pool.close()
            # simple progress bar
            while (True):
                remaining = rs._number_left
                print("finished:{0}/{1}".format(len(tasks) - remaining, len(tasks)),
                      end='\r')  # '\r' means remove the last line
                if (rs.ready()):
                    break
                time.sleep(0.5)

            pool.join()

        return None

    def show_psd(self):
        if os.path.exists(self.psd_figure_path):
            os.removedirs(self.psd_figure_path)
        os.makedirs(self.psd_figure_path)

        tele_df = pd.read_csv(self.tele_catalog)

        for tele_index, tele_row in tele_df.iterrows():
            tele_time = obspy.UTCDateTime(tele_row['time'])
            print('Plot %s...' % tele_time)

            tele_event_file = os.path.join(self.event_data_path_psd,
                                           str(tele_time),
                                           '.'.join([str(tele_time), self.sta, self.chn]))
            tb_b = obspy.UTCDateTime(tele_row['Tb_B_time'])
            tb_e = obspy.UTCDateTime(tele_row['Tb_E_time'])
            te_b = obspy.UTCDateTime(tele_row['Te_B_time'])
            te_e = obspy.UTCDateTime(tele_row['Te_E_time'])

            tele_tr = obspy.read(tele_event_file)[0]

            out_file = os.path.join(self.psd_figure_path, str(tele_time) + '.pdf')
            ax = plot_psd(tele_tr, tb_b, tb_e, te_b, te_e, self.time_segment, out_file=out_file)

        return None

    def show_bg_pir(self):
        if os.path.exists(self.background_pir_figure_path):
            os.removedirs(self.background_pir_figure_path)
        os.makedirs(self.background_pir_figure_path)

        tele_pir_file = os.path.join(
            self.tele_pir_path, self.sta + '.' + self.chn + '_PIR.csv')
        background_pir_associated_file = os.path.join(
            self.background_pir_associated_path,
            self.sta + '.' + self.chn + '_background_PIR_associated.csv')
        tele_pir_df = pd.read_csv(tele_pir_file)
        tele_background_pir_df = pd.read_csv(background_pir_associated_file)

        tele_pir_time = [obspy.UTCDateTime(t) for t in tele_pir_df['time'].values]

        for row_index, row in tele_background_pir_df.iterrows():
            tele_time = obspy.UTCDateTime(row['time'])

            if tele_time in tele_pir_time:
                pirs = re.split(' +', row['background_PIRs'][1:-1])
                pirs = [i for i in pirs if (i != '')]
                background_pir_list = np.array(pirs, dtype='float64')

                if len(background_pir_list) >= 30:
                    print (tele_time)
                    tele_pir = float(
                        tele_pir_df[tele_pir_df['time'] == row['time']]['PIR_log'].values)
                    out_file = os.path.join(self.background_pir_figure_path, str(tele_time) + '.pdf')
                    plot_bg_pir(tele_pir, background_pir_list, out_file=out_file)

    def cal_peak_dynamic_stress(self):
        if os.path.exists(self.peak_dynamic_stress_file):
            os.removedirs(self.peak_dynamic_stress_file)

        tele_df = pd.read_csv(self.tele_catalog)
        stress_array = []
        for tele_index, tele_row in tele_df.iterrows():
            tele_time = obspy.UTCDateTime(tele_row['time'])
            print('Plot %s...' % tele_time)

            tele_event_file = os.path.join(self.event_data_path_stress,
                                           str(tele_time),
                                           '.'.join([str(tele_time), self.sta, self.chn]))
            tele_event_tr = obspy.read(tele_event_file)[0]
            te_b = obspy.UTCDateTime(tele_row['Te_B_time'])
            te_e = obspy.UTCDateTime(tele_row['Te_E_time'])

            stress = peak_dynamic_stress(tele_event_tr,
                                         self.shear_wave_velocity,
                                         te_b, te_e,
                                         shear_modulus=self.shear_modulus)
            stress_array.append([str(tele_time), stress])

        stress_df = pd.DataFrame(data=stress_array,columns=['time', 'peak_dynamic_stress'])

        stress_df.to_csv(self.peak_dynamic_stress_file, index=False)

