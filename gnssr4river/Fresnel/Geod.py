# -*- coding: utf-8 -*-
"""
This code contains usefull functions for geodesic calculs.

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import numpy as np
import math as m

# Some parameters of WGS84 projection...
a = 6378137.0             # Earth radius (m)
f = 1./298.257223563      # flattening factor
e = np.sqrt(2*f-f**2)     # first eccentricity
b = a*(1-f)
e_pr = np.sqrt((a**2 - b**2) / b**2) # second excentricity

###############################################################################

# formules de passage du repère local au repère géocentrique

# (XL)   
# (YL) =  A . B 
# (ZL)

# Coordonnées du point M
# lambda1 = 
# phi1 = 

# A = np.array([       -np.sin(lambda1),               np.cos(lambda1),              0     ],
#                [-np.sin(phi1)*np.cos(lambda1),  -np.sin(phi1)*np.sin(lambda1), np.cos(phi1)],
#                [ np.cos(phi1)*np.cos(lambda1),   np.cos(phi1)*np.sin(lambda1), np.sin(phi1)])

# B = np.array([x - N1*np.cos(lambda1)*np.cos(phi1)],
#              [y - N1*np.sin(lambda1)*np.cos(phi1)],
#              [z - N1*(1-e**2)*np.sin(phi1)       ])

def GeotoCart(lon, lat, h):
    """
    This function transform geographical coordinates (WGS84) to cartesian coordinates (WGS84).
    
    Parameters
    ----------
    lat,lon: float
        location of latitude and longitudee in degrees 
    h: meters
        ellipsoïdal height

    Return
    ------
    X,Y,Z: float
        cartesian coordinates
        
    """
    # Put units to radians
    lon = m.radians(lon)
    lat = m.radians(lat)
    
    # Intermediate parameters
    w = 1 / np.sqrt(1 + e_pr**2*np.sin(lat)**2)
    N = a / w    
    
    # Coordinates
    X = (N + h)*np.cos(lon)*np.cos(lat)
    Y = (N + h)*np.sin(lon)*np.cos(lat)
    Z = (N*(1-e**2) + h)*np.sin(lat)
    
    return X,Y,Z

def CarttoGeo(X,Y,Z):
    """
    This function transform cartesian coordinates (WGS84) to geographical coordinates (WGS84).
    
    Parameters
    ----------
    X,Y,Z: float
        cartesian coordinates

    Return
    ------
    lat,lon: float
        location of latitude and longitudee in degrees 
    h: meters
        ellipsoïdal height
        
    """
    # Intermediate parameters
    f =  1 - np.sqrt(1-e**2)
    R = np.sqrt(X**2+Y**2+Z**2)
    mu = np.arctan((Z/np.sqrt(X**2+Y**2))*((1-f)+(e**2*a)/R))
    
    # Coordinates
    lon = np.arctan(Y/X)
    lat = np.arctan((Z*(1-f) + e**2*a*np.sin(mu)**3) / ((1-f)*(np.sqrt(X**2+Y**2) - e**2*a*np.cos(mu)**3)))
    h = np.sqrt(X**2+Y**2)*np.cos(lat) + Z*np.sin(lat) - a*np.sqrt(1 - e**2*np.sin(lat)**2)
    
    # Put units to degrees
    lon = m.degrees(lon)
    lat = m.degrees(lat)
    
    return lon, lat, h