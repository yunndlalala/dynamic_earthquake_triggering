#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/13
@file: example_main.py
"""
import json
import multiprocessing

from dyntripy.triggering import Triggering

import os
os.chdir(os.getcwd())

if __name__ == '__main__':
    # For Python 3.6 and 3.7, the start method of children processes should be "forkserver" or "spawn"
    multiprocessing.set_start_method("forkserver") # or 'spawn'

    # Load input file and define an instance of Triggering.
    hypers_f = open('input.json', 'r', encoding='utf-8')
    hypers = json.load(hypers_f)
    tri = Triggering(hypers)

    # Generate the power integral database with 2 processes.
    # tri.net_database(p=2)

    # Calculate the logarithmic ratios of remote earthquakes with 2 processes.
    # tri.net_ratio(p=2, bg=False)

    # Calculate the logarithmic ratios of background days with 2 processes.
    # tri.net_ratio(p=2, bg=True)

    # Compute the confidence levels with 2 processes.
    # tri.net_match(p=2)
    tri.net_cl(p=2)

    pass
