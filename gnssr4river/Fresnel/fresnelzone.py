# -*- coding: utf-8 -*-
"""
To compute and show the first Fresnel Zone.

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Usefull librairies
import numpy as np
import matplotlib.pyplot as plt
import math as m
from geod import *
from shapely.geometry.polygon import Polygon
import pyproj

# Usefull constants
c = 299792458                   # m.s-1 Speed of light
L1_GPS = 1575.42e6              # Hz L1 frequency for GPS
L2_GPS = 1227.60e6              # Hz L2 frequency for GPS
L1_Glo = 1602.0e6               # Hz L1 frequency for GLONASS
L2_Glo = 1246.0e6               # Hz 21 frequency for GLONASS
lambda_L1_GPS = (c/L1_GPS)      # m wavelenght for L1 GPS
lambda_L2_GPS = (c/L2_GPS)      # m wavelenght for L2 GPS
lambda_L1_Glo = (c/L1_Glo)      # m wavelenght for L1 Glo
lambda_L2_Glo = (c/L2_Glo)      # m wavelenght for L2 Glo

###############################################################################

# Calculation of the First Fresnel Zone

def firstFresnelZone(freq, h, elev):
    """
    This function gets the size and center of the First Fresnel Zone ellipse.
    (based on a code by Kristine Larson and Carolyn Roesler).
    
    Parameters
    ----------
    freq: float
        Frequence of L_band in Hz.
    h: float
        Hight of the receiver in meters.
    elev: float
        Satellite elevation angle in degrees.
    
    Return
    ------
    a: float
        Semi-major axis, aligned with the satellite azimuth (meters).
    b: float
        Semi-minor axis (meters).
    R: float 
        Locates the center of the ellispe on the satellite azimuth direction 
        and R meters away from the base of the Antenna.
    """
    # Valid frequencies
    lfreqgps = [L1_GPS, L2_GPS]
    lfreqglo = [L1_Glo, L2_Glo]

    # Check if frequency is valid
    if freq in lfreqgps and freq in lfreqglo:
        raise Exception("Wrong value for L_band frequency")  

    if elev>90:
        raise Exception("Wrong value for elevation, can't excede 90° !")  

    # Directly put elevation in radians
    elevR = np.radians(elev)
    
    if freq in lfreqgps:
        # delta = locus of points corresponding to a fixed delay;
        # typically the first Fresnel zone is is the
        # "zone for which the differential phase change across
        # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
        d = lambda_L1_GPS/2

        # from the appendix of Larson and Nievinski, 2013
        # semi-minor axis
        b = ((lambda_L1_GPS*h)/np.sin(elevR)) + (lambda_L1_GPS/(2*np.sin(elevR)))**2
        b = np.sqrt(b)
        # semi-majpr axis
        a = b/np.sin(elevR)

    elif freq in lfreqglo:
        # delta = locus of points corresponding to a fixed delay;
        # typically the first Fresnel zone is is the
        # "zone for which the differential phase change across
        # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
        d = lambda_L1_Glo/2

        # from the appendix of Larson and Nievinski, 2013
        # semi-minor axis
        b = ((lambda_L1_Glo*h)/np.sin(elevR)) + (lambda_L1_Glo/(2*np.sin(elevR)))**2
        b = np.sqrt(b)
        # semi-majpr axis
        a = b/np.sin(elevR)
        
    # determine distance to ellipse center in meters
    R = (h + d/np.sin(elevR)) / np.tan(elevR)

    return a, b, R

###############################################################################

def plotEllipse(a, b, R, lon, lat, h, azim):
    """
    Create an ellipse of a Fresnel zone.
    
    Parameters
    ----------
    a: float
        Semi-major axis, aligned with the satellite azimuth (meters).
    b: float
        Semi-minor axis (meters).
    R: float 
        Locates the center of the ellispe on the satellite azimuth direction 
        and R meters away from the base of the Antenna.
    lon,lat: float
        Position of the receiver in geographical coordinates (degrees).
    h: float
        Hight of the receiver in meters.
    azim: float
        Given azimut of ellipse in degrees.
    
    Return
    ------
    p: Polygon
        Polygon of the ellipse.
    area: float
        Area of the Polygon in square meter.
    """
    # Check input for the azimuth in case
    if azim > 360 or azim < 0:
        raise Exception("Wrong value of azimuth, should be between 0 and 360!")  
        
    # Set coordinates of the receiver to cartesians
    transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy = True)
    X = transformer.transform(lon, lat)[0]
    Y = transformer.transform(lon, lat)[1]
    
    # Set to radians
    azim = m.radians(azim)
    
    # Change angle to match orientation of Python
    angle = 2*np.pi - azim + np.pi/2
    
    # Coordinate of the center
    xR = X + R*np.cos(angle)                    
    yR = Y + R*np.sin(angle)
    
    t = np.linspace(0, 2*np.pi, 100)

    # Parametric equation of ellipse
    x = xR + a*np.cos(angle)*np.cos(t) - b*np.sin(angle)*np.sin(t)
    y = yR + a*np.sin(angle)*np.cos(t) + b*np.cos(angle)*np.sin(t)

    # Polygon of ellipse in epsg:3857
    q = Polygon(list(zip(x,y)))
    area = q.area
    
    # Changing back the coordinates to geographic
    lon = []
    lat = []
    for i in range(len(x)):
        lo = transformer.transform(x[i], y[i],direction='INVERSE')[0]
        la = transformer.transform(x[i], y[i],direction='INVERSE')[1]
        lon.append(lo)
        lat.append(la)
        
    # Polygon of ellipse in epsg:4326
    p = Polygon(list(zip(lon,lat)))

    return p, area

###############################################################################

def specularPoint(a, b, R, azim, color=None):
    """
    This function just return the center of an ellipse, i.e the reflection point.

    Parameters
    ----------
    a: float
       Semi-major axis, aligned with the satellite azimuth (meters).
    b: float
       Semi-minor axis (meters).
    R: float 
       Locates the center of the ellispe on the satellite azimuth direction 
       and R meters away from the base of the Antenna.
    azim: list 
       List of azimuths
    color: String (optional)
       Color of center points.

    Returns
    -------
    Plot of center of ellipses.

    """
    for angle in azim:
        
        angle = 2*np.pi - angle + np.pi/2
        xR = R*np.cos(angle)  # x-position of the center
        yR = R*np.sin(angle)  # y-position of the center

        if color != None:
            plt.axis('equal')
            plt.scatter(xR,yR, color=color)
        else:
            plt.axis('equal')
            plt.scatter(xR,yR)

    return

###############################################################################