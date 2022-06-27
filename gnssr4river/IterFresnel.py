# -*- coding: utf-8 -*-
"""
To iterate calculation of Fresnel Zone and show reflexion points for GNSS-R

@author: Lubin Roineau, ENSG-Geomatic (internship at UT-ITC Enschede), Aug 26, 2022
"""


# Usefull librairies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from FresnelZone import FresnelZone, PlotEllipse

# Usefull constants
c = 299792458  # m.s-1 Speed of light
L1_freq = 1575.42e6  # Hz L1 frequency
L2_freq = 1227.60e6  # Hz L2 frequency
lambda_L1 = c/L1_freq  # m wavelenght for L1
lambda_L2 = c/L2_freq  # m wavelenght for L2

