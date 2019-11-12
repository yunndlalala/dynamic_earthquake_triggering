'''
Written by Naidan YUN on 20180607.
'''
import numpy as np
import math
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")
def energy_attenuation():
    R=6371
    c_S = 4
    Q0 = 300
    yita=0.6
    FONTSIZE=15
    distance=[200,300,400,500,600,700,800,900,1000]
    f_list = np.arange(1, 50.1, 0.1)

    for d in distance:
        c1=1/(8*np.pi*(R**2))
        c2=1/math.sin(d/(2*R))
        EK_list=[]
        for f in f_list:
            Q=Q0*(f**yita)
            gama=(2 * np.pi * f) / (c_S*Q)
            EK=c1*c2*np.e**(-1*gama*d)
            # EK = np.e ** (-1 * gama * d)
            EK_list.append(EK)
        plt.plot(f_list, EK_list, '-')
        #print ('21')

    plt.legend(['distance=%dkm'%int(d) for d in distance])
    #plt.xlim([0,50])
    #plt.ylim([-0.01,0.9])

    tar_Energy_diff=10**(-9)
    plt.plot([0,50], [tar_Energy_diff, tar_Energy_diff], '--', color='k')

    plt.ylabel("E'/E",fontsize=FONTSIZE)
    plt.xlabel('f (Hz)',fontsize=FONTSIZE)
    plt.show()
    print ('Finish')
def travel_depth(d_in_degree):

    arrivals = model.get_ray_paths(source_depth_in_km=10,
                                   distance_in_degree=d_in_degree,
                                   phase_list=["S"])
    ax = arrivals.plot_rays(plot_type="cartesian")

if __name__ == '__main__':
    # distance = [80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    # d_in_degree=1000/111
    # travel_depth(d_in_degree)
    energy_attenuation()
    print ('Finish')

