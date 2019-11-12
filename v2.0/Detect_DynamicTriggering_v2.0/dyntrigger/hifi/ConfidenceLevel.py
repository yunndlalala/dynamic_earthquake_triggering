"""
@version:
author:yunnaidan
@time: 2018/11/23
@file: ConfidenceLevel.py
@function:
"""
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
#plt.rc('font',family='Times New Roman',size=16)
#plt.rc('font',family='Times New Roman')
plt.rc('font', family='Times New Roman', size=8)
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
import pandas as pd
from obspy.core import UTCDateTime
from datetime import timedelta
import seaborn as sbn
import statsmodels.api as sm
from scipy.stats import norm


def get_uniqe(list):
    out_list = [list[0]]
    out_list_index = [0]
    for i in range(1, len(list)):
        if list[i] not in out_list:
            out_list.append(list[i])
            out_list_index.append(i)
    return out_list


def read_associated_background_PIR(background_PIR_associated_file):
    with open(background_PIR_associated_file, 'r') as f:
        data = f.readlines()
    data_dict = {}
    for line in data:
        key = line.split(',')[0]
        value = [float('%.3f' % float(i)) for i in line[:-1].split(',')[1:]]
        data_dict[str(UTCDateTime(key))] = value
    return data_dict

def load_PIR(PIR_file, column_name):
    dataframe = pd.read_csv(PIR_file)
    events = get_uniqe(dataframe['event'].values)
    event_PIR_dict = {}
    for i in range(len(events)):
        event = events[i]
        PIR= dataframe[dataframe['event']== event][column_name].values[0]
        event_PIR_dict[event] = PIR
    return event_PIR_dict


def confidence_level_Gauss(
        background_PIR_associated_file,
        tele_PIR_file,
        out_file,
        PIR_columnname):

    tele_PIR_dict = load_PIR(tele_PIR_file, PIR_columnname)

    tele_background_PIR_dict = read_associated_background_PIR(
        background_PIR_associated_file)

    confidence_level_list = []

    for tele in tele_background_PIR_dict.keys():
        print('Calculate the confidence level of teleseism %s...' % str(tele))
        background_PIR_list = tele_background_PIR_dict[tele]
        x_norm, _,y_norm_cdf = sbn.distplot(
            background_PIR_list, fit=norm, kde=False)
        x = x_norm
        y = y_norm_cdf

        if tele in tele_PIR_dict.keys():
            tele_PIR = float(tele_PIR_dict[tele])
            for i in range(len(x)):
                if x[i] > tele_PIR:
                    left = x[i - 1]
                    right = x[i]
                    confidence_level_value = y[i - 1] + (((tele_PIR - left) / (
                        right - left)) * (y[i] - y[i - 1]))
                    break
                else:
                    confidence_level_value = 1
        else:
            confidence_level_value = None
        confidence_level_list.append([tele, confidence_level_value])

    CL_dataframe = pd.DataFrame(
        data=confidence_level_list, columns=[
            'tele', 'confidence_level_value'])
    CL_dataframe.to_csv(out_file, index=False)
    return CL_dataframe

def confidence_level(
        tele_backgroundHighPIR_sorted_file,
        tele_HighPIR_file,
        confidence_level_outfile):

    Tele_HighPIR_dict = load_PIR(tele_HighPIR_file, 'PIRatio')

    tele_backgroundTeleHighPIR_dict = read_associated_background_PIR(
        tele_backgroundHighPIR_sorted_file)
    tele_confidenceLevel_list = []
    for tele in tele_backgroundTeleHighPIR_dict.keys():
        #tele = '2001-07-07T09:38:43.520000Z'
        print('Calculate the confidence level of teleseism %s...' % str(tele))
        HighPID_list = tele_backgroundTeleHighPIR_dict[tele]
        kde = sm.nonparametric.KDEUnivariate(HighPID_list)
        kde.fit()
        x = kde.support
        y = 1 - kde.cdf
        #tele_key=tele[:-5] + 'Z'
        if tele in Tele_HighPIR_dict.keys():
            teleHighPID = np.log10(float(Tele_HighPIR_dict[tele]))
            for i in range(len(x)):
                if x[i] > teleHighPID:
                    left = x[i - 1]
                    right = x[i]
                    confidence_level_value = y[i - 1] + (((teleHighPID - left) / (
                        right - left)) * (y[i] - y[i - 1]))
                    break
                else:
                    confidence_level_value = 0
        else:
            confidence_level_value = None
        tele_confidenceLevel_list.append([tele, confidence_level_value])

    cdf_value_dataframe = pd.DataFrame(
        data=tele_confidenceLevel_list, columns=[
            'tele', 'confidence_level_value'])
    cdf_value_dataframe.to_csv(confidence_level_outfile, index=False)
    return cdf_value_dataframe

def view(
        bgPIR_list,
        tele_PIR,
        ax,
        hist=True,
        pdf=True,
        cdf=False,
        ylabel1=False,
        ylabel2=False,
        xlabel=False):

    tick_width = 0.5
    tick_length = 2

    if hist and pdf and (not cdf):
        x_norm, y_norm_pdf, y_norm_cdf = sbn.distplot(
            bgPIR_list, hist=True, fit=norm, kde=False, color="tab:gray", ax=ax, fit_kws={
                "color": "black", "label": 'PDF'})
        if xlabel:
            ax.set_xlabel('Logarithmic Ratio')
        if ylabel1:
            ax.set_ylabel('Probability')
        ylim = ax.get_ylim()
        ax.plot([np.mean(bgPIR_list), np.mean(bgPIR_list)],
                [ylim[0], ylim[1]], '--k')
        print('Mean logarithmic ratio is : %f' % np.mean(bgPIR_list))
        ax.plot([tele_PIR, tele_PIR], [
                ylim[0], ylim[1]], '--', color='tab:orange')

        ylim_new = ax.get_ylim()
        ax_pdf = ax.twinx()
        ax_pdf.plot(x_norm, y_norm_pdf, '-k')
        ax_pdf.set_ylim(ylim_new)
        ax_pdf.set_yticks([])

        ax.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)
        ax_pdf.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)
        mean = np.mean(bgPIR_list)
        var = (np.var(bgPIR_list)) ** 0.5
        # ax_pdf.set_xticks([])
        # ax.set_xticks([mean - 3 * var, mean - 2 * var, mean - var,
        #                mean, mean + var, mean + 2 * var, mean + 3 * var])
        # ax.set_xticklabels(
        #     ['0.001', '0.023', '0.159', '0.500', '0.841', '0.977', '0.999'])

        for i in range(len(x_norm)):
            if x_norm[i] > tele_PIR:
                left = x_norm[i - 1]
                right = x_norm[i]
                confidence_level_value = y_norm_cdf[i - 1] + (
                        ((tele_PIR - left) / (right - left)) * (y_norm_cdf[i] - y_norm_cdf[i - 1]))
                break
            else:
                confidence_level_value = 1
        print(confidence_level_value)
    if hist and pdf and cdf:
        x_norm, y_norm_pdf, y_norm_cdf = sbn.distplot(
            bgPIR_list, hist=True, fit=norm, kde=False, color="tab:gray", ax=ax, fit_kws={
                "color": "black", "label": 'PDF'})
        if xlabel:
            ax.set_xlabel('Logarithmic Ratio')
        if ylabel1:
            ax.set_ylabel('Normalized Number')
        ylim = ax.get_ylim()
        ax.plot([np.mean(bgPIR_list), np.mean(bgPIR_list)],
                [ylim[0], ylim[1]], '--k')
        #print ('Mean logarithmic ratio is : %f'%np.mean(bgHighPID_list))
        ax.plot([tele_PIR, tele_PIR], [
                ylim[0], ylim[1]], '--', color='tab:orange')

        ax_pdf = ax.twinx()
        # ax_pdf.set_xticks([])
        ax_pdf.plot(x_norm, y_norm_pdf, '-k')
        if ylabel2:
            ax_pdf.set_ylabel('Probability')

        ax_cdf = ax_pdf.twiny()
        ax_cdf.plot(y_norm_cdf, y_norm_pdf, '-k')

        ax.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)
        ax_pdf.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)

    if cdf and (not hist) and (not pdf):
        ax_hist = plt.figure().add_subplot()
        x_norm_cdf, y_norm_cdf = sbn.distplot(
            bgPIR_list, hist=False, fit=norm, kde=False, color="tab:gray", ax=ax_hist, fit_kws={
                "color": "black", "label": 'PDF'})
        ax.plot(x_norm_cdf, y_norm_cdf, color='black')

        ax.set_ylim([-0.1, 1.1])
        ax.plot([tele_PIR, tele_PIR],
                [-0.1, 1], '--', color='tab:orange')

        # Calculate and plot confidence level.
        for i in range(len(x_norm_cdf)):
            if x_norm_cdf[i] > tele_PIR:
                left = x_norm_cdf[i - 1]
                right = x_norm_cdf[i]
                confidence_level_value = y_norm_cdf[i - 1] + (
                        ((tele_PIR - left) / (right - left)) * (y_norm_cdf[i] - y_norm_cdf[i - 1]))
                break
            else:
                confidence_level_value = 1

        xlim = ax.get_xlim()
        ax.plot([xlim[0], xlim[1]], [confidence_level_value,
                                     confidence_level_value], '--', color='tab:orange')
        ax.set_xlim(xlim[0], xlim[1])
        ax.scatter(
            tele_PIR,
            confidence_level_value,
            c='tab:orange',
            edgecolors='tab:orange',
            marker='*',
            s=200)

        if xlabel:
            ax.set_xlabel('Logarithmic Ratio')
        if ylabel1:
            ax.set_ylabel('Confidence Level')

        ax.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)

    if (not hist) and pdf and cdf:
        x_norm_cdf, y_norm_cdf = sbn.distplot(
            bgPIR_list, hist=False, fit=norm, kde=False, color="tab:gray", ax=ax, fit_kws={
                "color": "black", "label": 'PDF'})

        if xlabel:
            ax.set_xlabel('Logarithmic Ratio')
        if ylabel1:
            ax.set_ylabel('Confidence Level')

        ax_cdf = ax.twinx()
        ax_cdf.plot(x_norm_cdf, y_norm_cdf, color='black')

        ylim = ax.get_ylim()
        ax.plot([tele_PIR, tele_PIR], [ylim[0], ylim[1]], '--r')

        ax.set_ylim([-0.1, 1.1])

        # Calculate and plot confidence level.
        for i in range(len(x_norm_cdf)):
            if x_norm_cdf[i] > tele_PIR:
                left = x_norm_cdf[i - 1]
                right = x_norm_cdf[i]
                confidence_level_value = y_norm_cdf[i - 1] + (
                        ((tele_PIR - left) / (right - left)) * (y_norm_cdf[i] - y_norm_cdf[i - 1]))
                break
            else:
                confidence_level_value = 1

        ax.scatter(
            tele_PIR,
            confidence_level_value,
            c='',
            edgecolors='tab:orange',
            marker='*',
            s=200)

        if xlabel:
            ax.set_xlabel('Logarithmic Ratio')
        if ylabel1:
            ax.set_ylabel('Confidence Level')

        ax.tick_params(
            which='both',
            direction='in',
            length=tick_length,
            width=tick_width)

    samples = len(bgPIR_list)
    return samples


if __name__ == '__main__':
    print('Finish!')
