#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/13
@file: merge_cl.py
"""
import os
import json
import platform

from dyntrigger.utils.preprocess_TeleCatalog import Bvel_Evel
from dyntrigger.hifi.HiFi import HiFi

if __name__ == '__main__':
    params_f = open('parameters.json', 'r', encoding='utf-8')
    params = json.load(params_f)
    stafile = os.path.join(ROOT_PATH,params['info']['sta_file'])
    sta_list = load_sta(stafile)

    # net_pi(
    #     sta_list,
    #     params['info']['chn_list'],
    #     params['PI']['time_segment'],
    #     params['PI']['f_win_list'],
    #     os.path.join(ROOT_PATH,params['PI']['Gf_info_file']),
    #     os.path.join(ROOT_PATH,params['PI']['raw_data_path']),
    #     os.path.join(ROOT_PATH,params['PI']['PIdatabase_path']),
    #     cores=15)
    #
    # net_tele_pir(
    #     sta_list,
    #     params['info']['chn_list'],
    #     os.path.join(ROOT_PATH,params['PI']['PIdatabase_path']),
    #     os.path.join(ROOT_PATH,params['tele_PIR']['tele_catalog']),
    #     os.path.join(ROOT_PATH,params['tele_PIR']['tele_PIR_path']),
    #     cores=1)
    #
    # Bvel_Evel(
    #     os.path.join(ROOT_PATH,params['tele_PIR']['tele_catalog']),
    #     params['background_PIR']['ref_lat'],
    #     params['background_PIR']['ref_lon'],
    #     params['background_PIR']['Tb'],
    #     params['background_PIR']['vel_B'],
    #     params['background_PIR']['vel_E'])
    # "tb": 18000,
    # "vel_b": 5,
    # "vel_e": 2,
    # "ref_lat": 35.81574,
    # "ref_lon": -117.59751,
    # catalog_during_days(
    #     os.path.join(ROOT_PATH,params['tele_PIR']['tele_catalog']),
    #     os.path.join(ROOT_PATH,params['background_PIR']['background_catalog_file']),
    #     dayWindow=params['background_PIR']['dayWindow'])
    # net_background_pir(
    #     sta_list,
    #     params['info']['chn_list'],
    #     os.path.join(ROOT_PATH,params['PI']['PIdatabase_path']),
    #     os.path.join(ROOT_PATH,params['background_PIR']['background_catalog_file']),
    #     os.path.join(ROOT_PATH,params['background_PIR']['background_PIR_path']),
    #     cores=1)
    # net_background_pir_asso(
    #     sta_list,
    #     params['info']['chn_list'],
    #     params['background_PIR']['dayWindow'],
    #     os.path.join(ROOT_PATH,params['tele_PIR']['tele_catalog']),
    #     os.path.join(ROOT_PATH,params['background_PIR']['background_PIR_path']),
    #     os.path.join(ROOT_PATH,params['background_PIR']['asso_bgPIR_path']),
    #     cores=1)

    net_cl(
        sta_list,
        params['info']['chn_list'],
        os.path.join(ROOT_PATH,params['tele_PIR']['tele_PIR_path']),
        os.path.join(ROOT_PATH,params['background_PIR']['asso_bgPIR_path']),
        os.path.join(ROOT_PATH,params['confidence_level']['CL_path']),
        cores=1)

    pass

