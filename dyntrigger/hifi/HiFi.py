#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: run_utils.py
"""
import os
import time
import multiprocessing

from dyntrigger.hifi import PIdatabase, HighPIRatio, AssociateBackgroundPIR, ConfidenceLevel
from dyntrigger.utils.preprocess_tele_catalog import catalog_during_days


class HiFi(object):

    def __init__(self, hypers):
        self.hypers = hypers
        self.sta_file = hypers['info']['sta_file']
        self.sta_list = self.load_sta()
        self.chn_list = hypers['info']['chn_list']
        self.time_segment = hypers['pi']['time_segment']
        self.f_win_list = hypers['pi']['f_win_list']
        self.day_window = hypers['background_PIR']['day_window']

        self.root_path = hypers['info']['root_path']
        self.gf_info = os.path.join(self.root_path, hypers['pi']['gf_info'])
        self.raw_data_path = os.path.join(self.root_path, hypers['pi']['raw_data_path'])
        self.pi_database_path = os.path.join(self.root_path, hypers['pi']['pi_database_path'])
        self.tele_catalog = os.path.join(self.root_pat, hypers['tele_pir']['tele_catalog'])
        self.tele_pir_path = os.path.join(self.root_path, hypers['tele_pir']['tele_pir_path'])
        self.background_catalog = os.path.join(self.root_path, hypers['background_pir']['background_catalog'])
        self.background_pir_path = os.path.join(self.root_path, hypers['background_pir']['background_pir_path'])
        self.background_pir_associated_path = os.path.join(self.root_path,
                                                           hypers['background_pir']['background_pir_associated_path'])
        self.cl_path = os.path.join(self.root_path, hypers['confidence_level']['cl_path'])

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

    def net_pi(self, cores):
        if os.path.exists(self.pi_database_path):
            os.removedirs(self.pi_database_path)
        os.makedirs(self.pi_database_path)

        for sta in self.sta_list:
            print(sta)
            for chn in self.chn_list:
                print(' ' + chn)
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
                print(sta)
                for chn in self.chn_list:
                    print(' ' + chn)
                    pi_database_folder = os.path.join(self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(self.tele_pir_path, sta + '.' + chn + '_PIR.csv')
                    HighPIRatio.run_pir(
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
                        (self.tele_catalog,
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

    def net_background_pir(self,cores=1):

        catalog_during_days(self.tele_catalog, self.background_catalog, self.day_window)

        if os.path.exists(self.background_pir_path):
            os.removedirs(self.background_pir_path)
        os.makedirs(self.background_pir_path)

        if cores == 1:
            for sta in self.sta_list:
                print(sta)
                for chn in self.chn_list:
                    print(' ' + chn)
                    pi_database_folder = os.path.join(
                        self.pi_database_path, sta + '.' + chn)
                    out_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    HighPIRatio.run_pir(
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
                        (self.background_catalog,
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
                print(sta)
                for chn in self.chn_list:
                    print(' ' + chn)
                    background_pir_file = os.path.join(
                        self.background_pir_path, sta + '.' + chn + '_background_PIR.csv')
                    out_file = os.path.join(
                        self.background_pir_associated_path, sta + '.' + chn + '_background_PIR_associated.csv')
                    AssociateBackgroundPIR.run_associate(
                        self.tele_catalog, background_pir_file, 'PIRatio_log', self.day_window, out_file)
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
                         'PIRatio_log',
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
                print(sta)
                for chn in self.chn_list:
                    print(' ' + chn)
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
                        'PIRatio_log')
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
                        'PIRatio_log'))
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
