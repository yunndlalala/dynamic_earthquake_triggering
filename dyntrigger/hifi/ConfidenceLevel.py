"""
author:yunnaidan
@time: 2018/11/23
@file: ConfidenceLevel.py
"""
import numpy as np
import pandas as pd
from scipy import stats
import re


def cl41event(row, PIR_columnname, tele_PIR_df):
    if row['time'] in tele_PIR_df['time'].values:
        PIRs = re.split(' +', row['background_PIRs'][1:-1])
        PIRs = [i for i in PIRs if (i != '')]
        background_PIR_list = np.array(PIRs, dtype='float64')
        if len(background_PIR_list) < 30:
            confidence_level_value = 'missing background PIR'
        else:
            tele_PIR = float(
                tele_PIR_df[tele_PIR_df['time'] == row['time']][PIR_columnname].values)

            mean, std = stats.norm.fit(background_PIR_list)
            background_norm = stats.norm(loc=mean, scale=std)

            confidence_level_value = background_norm.cdf(tele_PIR)
    else:
        confidence_level_value = 'missing teleseismic PIR'

    return confidence_level_value


def run_cl(
        background_PIR_associated_file,
        tele_PIR_file,
        out_file,
        PIR_columnname):

    tele_PIR_df = pd.read_csv(tele_PIR_file)

    tele_background_PIR_df = pd.read_csv(background_PIR_associated_file)

    confidence_level_list = []

    for row_index, row in tele_background_PIR_df.iterrows():
        print ('Calculate cl of %s' % row['time'])
        confidence_level_value = cl41event(row, PIR_columnname, tele_PIR_df)
        confidence_level_list.append([row['time'], confidence_level_value])

    CL_dataframe = pd.DataFrame(
        data=confidence_level_list, columns=[
            'time', 'confidence_level_value'])
    CL_dataframe.to_csv(out_file, index=False)

    return CL_dataframe
