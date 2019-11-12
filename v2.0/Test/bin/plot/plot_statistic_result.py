"""
@version:
author:yunnaidan
@time: 2019/03/07
@file: Figure5.py
@function:
"""
import os
import matplotlib.pyplot as plt
from dyntrigger.hifi.plot_utils.visual_hifiResult import background_distribution_plot
plt.rc('font', family='Times New Roman', size=8)
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
# plot paras
tick_width = 0.5
tick_length = 2
# Parameters
sta='NC.GDXB'
chn='HHZ'
tele = '2010-04-04T22:40:42.36Z'
tele_PIR_file = '../../result/tele_PIR/NC.GDXB.HHZ_PIR.csv'
background_PIR_associated_path = '../../result/background_PIR_associated/NC.GDXB.HHZ_background_PIR_associated.csv'
kwags = {'hist': True, 'pdf': True, 'cdf': False,'ylabel1':True,'ylabel2':True,'xlabel':True}

# Plot
fig = plt.figure()

ax = fig.add_subplot(111)
background_distribution_plot(
    tele_PIR_file,
    background_PIR_associated_path,
    tele,
    ax, **kwags)

ax.tick_params(
    which='both',
    direction='in',
    length=tick_length,
    width=tick_width)

# out_figfile = '/Users/yunnaidan/Project/Dynamic_Triggering/Documents/MINE/Paper/Figures_Tables/version3/Figure5'
# plt.savefig(os.path.join(out_figfile,'figure5_sup.pdf'))
plt.show()
print('Finish')
