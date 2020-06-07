#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: yunnaidan
@time: 2019/11/13
@file: example.py
"""
import json
from dyntripy.triggering import Triggering

# Load input file and define an instance of Triggering.
hypers_f = open('input.json', 'r', encoding='utf-8')
hypers = json.load(hypers_f)
tri = Triggering(hypers)

# Generate the power integral database with two processes.
tri.net_database(p=2)

# Calculate the logarithmic ratios of remote earthquakes with two processes.
tri.net_ratio(p=2, bg=False)

# Calculate the logarithmic ratios of background days with two processes.
tri.net_ratio(p=2, bg=True)

# Compute the confidence levels with two processes.
tri.net_cl(p=2)