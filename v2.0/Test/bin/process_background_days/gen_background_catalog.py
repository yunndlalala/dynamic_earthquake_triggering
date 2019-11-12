"""
@version:
author:yunnaidan
@time: 2018/10/08
@file: gen_background_catalog.py
@function:
"""
from dyntrigger.utils.pyTOOL.preprocess_TeleCatalog import catalog_during_days

def main():
    # Catalog of the telesismic events
    teleseismic_catalog='../../data/teleseismic_catalog.csv'
    # Catalog of the background days
    out_file='../../data/background_catalog.csv'
    # The number of background days for one teleseismic event
    dayWindow=4 # 120

    catalog_during_days(teleseismic_catalog, out_file,dayWindow=dayWindow)

if __name__ == "__main__":
    main()
    print ('Finish!')
