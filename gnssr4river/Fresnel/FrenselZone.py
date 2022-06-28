# -*- coding: utf-8 -*-
"""
To compute and show the first Fresnel Zone

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Usefull librairies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Usefull constants
c = 299792458  # m.s-1 Speed of light
L1_freq = 1575.42e6  # Hz L1 frequency
L2_freq = 1227.60e6  # Hz L2 frequency
lambda_L1 = c/L1_freq  # m wavelenght for L1
lambda_L2 = c/L2_freq  # m wavelenght for L2


# Calculation of the First Fresnel Zone

def FresnelZone(freq, h, elev, azim):
    """
    This function gets the size and center of the First Fresnel Zone ellipse
    at the selected  L-band frequency freq (Hz)
    for an Antenna height h (meters) above the reflecting surface
    for a satellite elevation angle elev (degrees)
    and azimuth direction azim (degrees)

    Output are the parameters of the ellipse [a, b, R ] (meters) where:
    a is the semi-major axis, aligned with the satellite azimuth 
    b is the semi-minor axis
    R locates the center of the ellispe 
    on the satellite azimuth direction (azim)
    and R meters away from the base of the Antenna.

    The ellipse is located on a flat horizontal surface h meters below
    the receiver.                  

    (based on a code by Kristine Larson and Carolyn Roesler)
    """

    lfreq = [L1_freq, L2_freq]

    if freq not in lfreq:
        return "Wrong value for L_band frequency"

    # delta = locus of points corresponding to a fixed delay;
    # typically the first Fresnel zone is is the
    # "zone for which the differential phase change across
    # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
    d = c/freq/2

    # semi-major and semi-minor dimension
    # from the appendix of Larson and Nievinski, 2013
    sin_elev_deg = np.sin(elev) * 180/np.pi
    b = (2*d*h/sin_elev_deg) + (d/sin_elev_deg)**2
    if b < 0:
        return "Invalid value for b"
    b = np.sqrt(b)
    a = b/sin_elev_deg

    # determine distance to ellipse center
    R = (h + d/(np.sin(elev) * 180/np.pi)) / (np.tan(elev) * 180/np.pi)

    return a, b, R


def plotEllipse(a, b, R, azim):
    """
    Plots an ellipse
    
    a is the semi-major axis, aligned with the satellite azimuth 
    b is the semi-minor axis
    R locates the center of the ellispe (distance in meters from receiver)
    on the satellite azimuth direction (azim)

    We assume the receiver is located within a grid at coordinates (0,0)
    """

    xR = R*np.sin(azim)  # x-position of the center
    yR = R*np.cos(azim)  # y-position of the center
    rot = azim  # rotation angle

    t = np.linspace(0, 2*np.pi, 100)
    Ell = np.array([a*np.cos(t), b*np.sin(t)])
    # Rotation matrix
    R_rot = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])

    Ell_rot = np.zeros((2, Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:, i] = np.dot(R_rot, Ell[:, i])

    plt.plot(xR+Ell_rot[0, :], yR+Ell_rot[1, :])  # rotated ellipse
    plt.show()
