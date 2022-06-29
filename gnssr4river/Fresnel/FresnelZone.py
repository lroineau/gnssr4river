# -*- coding: utf-8 -*-
"""
To compute and show the first Fresnel Zone

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Usefull librairies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import math as m

# Usefull constants
c = 299792458  # m.s-1 Speed of light
L1_freq = 1575.42e6  # Hz L1 frequency
L2_freq = 1227.60e6  # Hz L2 frequency
lambda_L1 = (c/L1_freq)  # m wavelenght for L1
lambda_L2 = (c/L2_freq)  # m wavelenght for L2


# Calculation of the First Fresnel Zone

def FirstFresnelZone(freq, h, elev):
    """
    This function gets the size and center of the First Fresnel Zone ellipse
    at the selected  L-band frequency freq (Hz)
    for an Antenna height h (meters) above the reflecting surface
    for a satellite elevation angle elev (degrees)

    Output are the parameters of the ellipse [a, b, R ] (meters) where:
    b is the semi-major axis, aligned with the satellite azimuth 
    a is the semi-minor axis
    R locates the center of the ellispe 
    on the satellite azimuth direction
    and R meters away from the base of the Antenna.

    The ellipse is located on a flat horizontal surface h meters below
    the receiver.                  

    (based on a code by Kristine Larson and Carolyn Roesler)
    """

    lfreq = [L1_freq, L2_freq]

    if freq not in lfreq:
        raise Exception("Wrong value for L_band frequency")  

    if elev>50:
        raise Exception("Wrong value for elevation, can't excede 50° !")  


    # delta = locus of points corresponding to a fixed delay;
    # typically the first Fresnel zone is is the
    # "zone for which the differential phase change across
    # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
    d = lambda_L1/2

    # from the appendix of Larson and Nievinski, 2013
    # semi-minor axis
    b = ((lambda_L1*h)/np.sin(m.radians(elev))) + (lambda_L1/(2*np.sin(m.radians(elev))))**2
    b = np.sqrt(b)
    # semi-majpr axis
    a = b/np.sin(m.radians(elev))

    # determine distance to ellipse center
    R = (h + d/(np.sin(elev) * np.pi/180)) / (np.tan(elev) * np.pi/180)

    return a, b, R


def plotEllipse(a, b, R, azim):
    """
    Plots an ellipse
    
    a is the semi-major axis, aligned with the satellite azimuth 
    b is the semi-minor axis
    R locates the center of the ellispe (distance in meters from receiver)
    on the satellite azimuth direction (azim)

    We assume the receiver is located at coordinates (0,0), y axis for North and x axis for East
    """

    azim = azim *np.pi/180
    xR = R*np.sin(azim)  # x-position of the center
    yR = R*np.cos(azim)  # y-position of the center
    theta = (90*np.pi/180)+azim  # rotation angle

    t = np.linspace(0, 2*np.pi, 100)

    x = xR + a*np.cos(theta)*np.cos(t) - b*np.sin(theta)*np.sin(t)
    y = yR + a*np.sin(theta)*np.cos(t) + b*np.cos(theta)*np.sin(t)
    # plt.axis('equal')
    plt.plot(x, y)  # rotated ellipse


def Center(a, b, R):
    azim = np.linspace(0, 2*np.pi, 50)    
    for angle in azim:
        xR = R*np.sin(angle)  # x-position of the center
        yR = R*np.cos(angle)  # y-position of the 
        plt.scatter(xR,yR)
    plt.show()