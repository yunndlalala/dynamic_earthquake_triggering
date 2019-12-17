#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: gf.py
"""
import numpy as np
import matplotlib.pyplot as plt

from dyntrigger.hifi.PIdatabase import gf


def plot_Gf():
    sensivity = 1
    zeros = [0, 0, 0, 0]
    poles = [-5.0265 + 3.76987j,
             -5.0265 - 3.76987j,
             -0.5969,
             -0.5969,
             -276.46,
             -276.46,
             -48.0915 + 116.097j,
             -48.0915 - 116.097j,
             -116.101 + 48.0832j,
             -116.101 - 48.0832j]
    normalizing = 1.936754 * (10 ** 13)
    f_list = np.arange(0.001, 1000, 0.001)
    Gf_list = []
    for f in f_list:
        Gf_list.append(gf(sensivity, normalizing, zeros, poles, f))
    print(gf(sensivity, normalizing, zeros, poles, 1))
    plt.plot(f_list, Gf_list)
    plt.plot([0.001, 100], [1, 1])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
