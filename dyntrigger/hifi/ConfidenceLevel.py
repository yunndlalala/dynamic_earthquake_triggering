"""
author:yunnaidan
@time: 2018/11/23
@file: ConfidenceLevel.py
"""
import numpy as np
import pandas as pd
from scipy import stats
import re
import logging

logging.basicConfig(
    filename='log.txt',
    format='%(asctime)s-%(name)s-%(levelname)s-%(module)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S %p',
    level=logging.INFO,
    filemode='w')


def cl41event(row, pir_columnname, tele_pir_df):
    if row['time'] in tele_pir_df['time'].values:
        pirs = re.split(' +', row['background_PIRs'][1:-1])
        pirs = [i for i in pirs if (i != '')]
        background_pir_list = np.array(pirs, dtype='float64')
        if len(background_pir_list) < 30:
            confidence_level_value = 'Insufficient background PIR'
        else:
            tele_pir = float(
                tele_pir_df[tele_pir_df['time'] == row['time']][pir_columnname].values)

            mean, std = stats.norm.fit(background_pir_list)
            background_norm = stats.norm(loc=mean, scale=std)

            confidence_level_value = background_norm.cdf(tele_pir)
    else:
        confidence_level_value = 'No teleseismic PIR'

    return confidence_level_value


def run_cl(
        background_pir_associated_file,
        tele_pir_file,
        out_file,
        pir_column_name):

    tele_pir_df = pd.read_csv(tele_pir_file)

    tele_background_pir_df = pd.read_csv(background_pir_associated_file)

    confidence_level_list = []

    for row_index, row in tele_background_pir_df.iterrows():
        print ('Calculate cl of %s' % row['time'])
        confidence_level_value = cl41event(row, pir_column_name, tele_pir_df)
        confidence_level_list.append([row['time'], confidence_level_value])

    cl_dataframe = pd.DataFrame(
        data=confidence_level_list, columns=[
            'time', 'confidence_level_value'])
    cl_dataframe.to_csv(out_file, index=False)

    return None
