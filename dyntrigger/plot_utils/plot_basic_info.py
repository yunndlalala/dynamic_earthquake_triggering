#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/12
@file: gf.py
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', size=10)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

from dyntrigger.utils.basic_utils import load_gf, gf


def plot_gf(
        sac_file, gf_info_file,
        frequencies=np.arange(0.01, 20, 0.0001), ax=None, **plt_paras):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    type, sensitivity, normalizing, zeros, poles = load_gf(sac_file, gf_info_file)
    gf_list = np.array(
        [gf(sensitivity, normalizing, zeros, poles, f) for f in frequencies])

    ax.plot(frequencies, gf_list, **plt_paras)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Frequency ' + '$\mathregular{(Hz)}$')
    ax.set_ylabel('Count/' + type)

    return ax


