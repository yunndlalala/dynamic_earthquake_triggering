"""
@version:
author:yunnaidan
@time: 2018/10/26
@file: time_window.py
@function:
"""
from dyntrigger.utils.pyTOOL.preprocess_TeleCatalog import Bvel_Evel



def main():
    # Time window before the arrival of the teleseismic event
    Tb=3600*5 # Seconds
    # Start of the time window after the arrival
    vel_B = 5 # km/s
    # End of the time window after the arrival
    vel_E = 2 # km/s
    # latitude of station
    sta_lat=38.8080
    # longitude of station
    sta_lon=-122.7953
    # File of catalog
    catalog_file = '../data/teleseismic_catalog.csv'

    Bvel_Evel(catalog_file, sta_lat,sta_lon,Tb,vel_B, vel_E)
if __name__ == '__main__':
    main()
    print ('Finish!')
