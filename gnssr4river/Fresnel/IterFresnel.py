# -*- coding: utf-8 -*-
"""
To iterate calculation of Fresnel Zone and show reflexion points for GNSS-R

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""


# Usefull librairies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from FresnelZone import FirstFresnelZone, plotEllipse
from GetOrbits import *

# Usefull constants
c = 299792458  # m.s-1 Speed of light
L1_freq = 1575.42e6  # Hz L1 frequency
L2_freq = 1227.60e6  # Hz L2 frequency
lambda_L1 = c/L1_freq  # m wavelenght for L1
lambda_L2 = c/L2_freq  # m wavelenght for L2

elev = [5,10,15]
azim = np.linspace(0, 360, 1000)
a5,b5,R5 = FirstFresnelZone(L1_freq, 2, 5)
a10,b10,R10 = FirstFresnelZone(L1_freq, 2, 10)
a15,b15,R15 = FirstFresnelZone(L1_freq, 2, 15)


for angle in azim:
    plotEllipse(a5,b5,R5,angle)
    plotEllipse(a10,b10,R10,angle)
    plotEllipse(a15,b15,R15,angle)

    
plt.xlabel("E (m)")
plt.ylabel("N (m)")
plt.title("Fresnel Surface")
plt.legend()
plt.show()        