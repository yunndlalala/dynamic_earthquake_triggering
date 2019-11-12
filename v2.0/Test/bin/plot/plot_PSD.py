import os
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman', size=8)
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
import obspy
from obspy import UTCDateTime
from pylab import *
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
from dyntrigger.hifi.plot_utils.visual_hifiResult import powerSpectralDensity_plot

# plot para
tick_width = 0.5
tick_length = 2

tele = '2010-04-04T22:40:42.36Z'
raw_data_path='../../data/raw/2010/20100404'

Tb_B = UTCDateTime('2010-04-04T17:42:55.420781Z')
Tb_E = UTCDateTime('2010-04-04T22:42:55.420781Z')
Te_B = UTCDateTime('2010-04-04T22:44:07.140000Z')
Te_E = UTCDateTime('2010-04-04T22:49:14.310000Z')

sta='NC.GDXB'
chn='HHZ'

sensivity = 1.006000e+09
zeros = [0,
         0,
         complex(-4.618140e+02, -4.290780e+02),
         complex(-4.618140e+02, 4.290780e+02),
         complex(-1.946970e+02, 0),
         complex(-1.514870e+01, 0)]
poles = [complex(-3.681800e-02, -3.692300e-02),
         complex(-3.681800e-02, 3.692300e-02),
         complex(-1.023970e+04, -2.725020e+03),
         complex(-1.023970e+04, 2.725020e+03),
         complex(-9.512740e+03, -1.146990e+04),
         complex(-9.512740e+03, 1.146990e+04),
         complex(-4.545250e+02, 0.000000e+00),
         complex(-3.981840e+02, 0.000000e+00),
         complex(-9.294710e+01, -3.889920e+02),
         complex(-9.294710e+01, 3.889920e+02),
         complex(-1.545030e+01, 0.000000e+00)]
normalizing = 9.482280e+18
Gf_parameters = [sensivity, normalizing, zeros, poles]

time_segment=30


####################
fig = plt.figure()

ax = fig.add_subplot(111)
raw_data_file=os.path.join(raw_data_path,sta+'.'+chn)
day_st = obspy.read(raw_data_file)
day_tr = day_st[0]
powerSpectralDensity_plot(
    day_tr,
    Tb_B,
    Tb_E,
    Te_B,
    Te_E,
    Gf_parameters,
    time_segment,
    ax)

ax.set_ylabel(
    'PSD ' +
    '$\mathregular{((nm/s)^2/(Hz\cdot s))}$')
ax.set_xlabel('Frequency (Hz)')
ax.plot([25, 25], [10**9,10**9], '--k')
ax.plot([35, 35], [10**9,10**9], '--k')

ax.tick_params(
which='both',
direction='in',
length=tick_length,
width=tick_width)
ax.set_xlim([1, 40])
ax.set_ylim([10,10**11])





# outpath = '/Users/yunnaidan/Project/Dynamic_Triggering/Documents/MINE/Paper/Figures_Tables/version_GRL/Figure1'
# plt.savefig(os.path.join(outpath, 'figure1.pdf'))
fig.show()
print('Finish!')
