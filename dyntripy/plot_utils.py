#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2020/08/09
@file: plot_utils.py
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def distribution(
        background_pir_list,
        tele_pir,
        out_file,
        ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    mean, std = stats.norm.fit(background_pir_list)
    background_norm = stats.norm(loc=mean, scale=std)
    confidence_level_value = background_norm.cdf(tele_pir)
    x = np.linspace(
        background_norm.ppf(0.001),
        background_norm.ppf(0.999),
        100)
    y_pdf = background_norm.pdf(x)

    ax.hist(
        background_pir_list,
        density=True,
        alpha=0.5,
        color='gray')
    ax.plot(x, y_pdf, '-', color='gray')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.plot([mean, mean], [ymin, ymax], '-k')
    ax.plot([tele_pir, tele_pir], [ymin, ymax], '--', color='k')
    ax.text(xmin + 0.8 * (xmax - xmin),
            ymin + 0.9 * (ymax - ymin),
            'CL = %.3f' % confidence_level_value)
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel('Logarithmic Ratio')
    ax.set_ylabel('Probability Density')

    ax_top = ax.twiny()
    ax_top.set_xticks([mean - 3 * std,
                       mean - 2 * std,
                       mean - std,
                       mean,
                       mean + std,
                       mean + 2 * std,
                       mean + 3 * std])
    ax_top.set_xticklabels(
        ['0.001', '0.023', '0.159', '0.500', '0.841', '0.977', '0.999'])
    ax_top.set_xlabel('Confidence Level')
    xlim = ax.get_xlim()
    ax_top.set_xlim(xlim)

    plt.savefig(out_file)

    plt.close()

    return None
