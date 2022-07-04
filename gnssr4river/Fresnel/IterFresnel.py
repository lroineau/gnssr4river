# -*- coding: utf-8 -*-
"""
To iterate calculation of Fresnel Zone and show reflexion points for GNSS-R

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede, NL), Aug 26, 2022
"""


# Usefull librairies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from FresnelZone import FirstFresnelZone, plotEllipse
from GetOrbits import *
from PlotFresnel import *

# Usefull constants
c = 299792458             # m.s-1 Speed of light
L1_freq = 1575.42e6       # Hz L1 frequency
L2_freq = 1227.60e6       # Hz L2 frequency
lambda_L1 = c/L1_freq     # m wavelenght for L1
lambda_L2 = c/L2_freq     # m wavelenght for L2

def IterFresnel(freq, h, elev_list, azim_min, azim_max, nb=50):
    """
    Iterate the FirstFresnelZone calcule.
    
    Parameters
    ----------
    freq: float
        frequence of l_band
    h: float
        hight of the receiver
    elev: list
        list of elevation angle in degrees (should not exceed 4 values for vizualisation)
    azim_min, azim_max: float
        minimum and maximum azimut for calcule  
    nb: int
        number of ellipse to plot
    
    Return
    ------
    A plot with all Fresnel zones
    """
    # Set azimut and color
    azim = np.linspace(azim_min, azim_max, nb)
    color = ['yellow','blue','red']
    
    # Check if elevation is a list or a single value
    if isinstance(elev_list, list):
        print('More than one elevation given, processing...')
        # Calcule First Fresnel Surface for each elevation angles
        c = 0
        for elev in elev_list:
            for angle in azim:
                a,b,R = FirstFresnelZone(freq, h, elev)
                plotEllipse(a,b,R,angle, color[c])
            c+=1
                    
    else:
        elev = elev_list
        print('only one elevation given, processing...')
        # Calcule First Fresnel Surface for each elevation angle
        for angle in azim:
            a,b,R = FirstFresnelZone(freq, h, elev)
            plotEllipse(a,b,R,angle, 'green')
                
    return        
