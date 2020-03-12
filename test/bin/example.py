#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/13
@file: example.py
"""
import os
import json

from dyntrigger.hifi.HiFi import HiFi
from dyntrigger.utils.preprocess_tele_catalog import begin_v_end_v


hypers_f = open('parameters.json', 'r', encoding='utf-8')
hypers = json.load(hypers_f)

detector = HiFi(hypers)

my_cores = 2
detector.net_pi(cores=my_cores)

Tb = 18000
vel_B = 5.0
vel_E = 2.0
ref_lat = 35.81574
ref_lon = -117.59751
begin_v_end_v(os.path.join(hypers['info']['root_path'],hypers['tele_pir']['tele_catalog']),
          ref_lat, ref_lon, Tb, vel_B, vel_E)

detector.net_tele_pir(cores=my_cores)

detector.gen_background_catalog()
detector.net_background_pir(cores=my_cores)

detector.net_background_pir_asso(cores=my_cores)

detector.net_cl(cores=my_cores)

detector.show_psd()
detector.show_bg_pir()

detector.cal_peak_dynamic_stress()

pass
