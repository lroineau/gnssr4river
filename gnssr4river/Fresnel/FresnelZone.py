# -*- coding: utf-8 -*-
"""
To compute and show the first Fresnel Zone

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Usefull librairies
import numpy as np
import matplotlib.pyplot as plt
import math as m
from Geod import CarttoGeo, GeotoCart

# Usefull constants
c = 299792458               # m.s-1 Speed of light
L1_GPS = 1575.42e6         # Hz L1 frequency for GPS
L2_GPS = 1227.60e6         # Hz L2 frequency for GPS
L1_Glo = 1575.42e6         # Hz L1 frequency for GPS
L2_Glo = 1227.60e6         # Hz L2 frequency for GPS
lambda_L1_GPS = (c/L1_GPS)     # m wavelenght for L1 GPS
lambda_L2_GPS = (c/L2_GPS)     # m wavelenght for L2 GPS
lambda_L1_Glo = (c/L1_Glo)     # m wavelenght for L1 Glonass
lambda_L2_Glo = (c/L2_Glo)     # m wavelenght for L2 Glonass



###############################################################################

# Calculation of the First Fresnel Zone

def FirstFresnelZone(freq, h, elev):
    """
    This function gets the size and center of the First Fresnel Zone ellipse.
    (based on a code by Kristine Larson and Carolyn Roesler)
    
    Parameters
    ----------
    freq: float
        frequence of l_band in Hz
    h: float
        hight of the receiver in meters
    elev: float
        satellite elevation angle in degrees
    
    Return
    ------
    a: float
        semi-major axis, aligned with the satellite azimuth (meters)
    b: float
        semi-minor axis (meters)
    R: float 
        locates the center of the ellispe on the satellite azimuth direction 
        and R meters away from the base of the Antenna.
    """

    lfreqgps = [L1_GPS, L2_GPS]
    lfreqglo = [L1_Glo, L2_Glo]

    if freq not in lfreqgps or freq not in lfreqglo:
        raise Exception("Wrong value for L_band frequency")  

    if elev>90:
        raise Exception("Wrong value for elevation, can't excede 90° !")  

    if freq in lfreqgps:
        # delta = locus of points corresponding to a fixed delay;
        # typically the first Fresnel zone is is the
        # "zone for which the differential phase change across
        # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
        d = lambda_L1_GPS/2

        # from the appendix of Larson and Nievinski, 2013
        # semi-minor axis
        b = ((lambda_L1_GPS*h)/np.sin(m.radians(elev))) + (lambda_L1_GPS/(2*np.sin(m.radians(elev))))**2
        b = np.sqrt(b)
        # semi-majpr axis
        a = b/np.sin(m.radians(elev))

    elif freq in lfreqglo:
        # delta = locus of points corresponding to a fixed delay;
        # typically the first Fresnel zone is is the
        # "zone for which the differential phase change across
        # the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
        d = lambda_L1_GPS/2

        # from the appendix of Larson and Nievinski, 2013
        # semi-minor axis
        b = ((lambda_L1_GPS*h)/np.sin(m.radians(elev))) + (lambda_L1_GPS/(2*np.sin(m.radians(elev))))**2
        b = np.sqrt(b)
        # semi-majpr axis
        a = b/np.sin(m.radians(elev))
        
    # determine distance to ellipse center in meters
    R = (h + d/np.sin(m.radians(elev))) / np.tan(m.radians(elev))

    return a, b, R

###############################################################################

def plotEllipse(a, b, R, lon, lat, h, azim, color):
    """
    Plot an ellipse
    
    Parameters
    ----------
    a: float
        semi-major axis, aligned with the satellite azimuth (meters)
    b: float
        semi-minor axis (meters)
    R: float 
        locates the center of the ellispe on the satellite azimuth direction 
        and R meters away from the base of the Antenna.
    lon,lat: float
        position of the receiver in geographical coordinates (degrees)
    h: float
        hight of the receiver in meters
    azim: float
        given azimut of ellipse in degrees
    color: String
        given color for the plot
    
    Return
    ------
    A plot of the ellipse
    """
    if azim > 360 or azim < 0:
        raise Exception("Wrong value azim, should be between 0 and 360!")  
        
    # Set coordinates of the receiver to cartesians
    X,Y,Z = GeotoCart(lon,lat,h)

    azim = azim *np.pi/180
    #xR = R*np.sin(azim)  # x-position of the center
    #yR = R*np.cos(azim)  # y-position of the center
    theta = 0 
        
    if 0<=azim<=np.pi/2:
        theta = (90*np.pi/180) - azim  # rotation angle
#        print('theta angle;', theta)
        xR = X + R*np.sin(theta)  # x-position of the center
        yR = Y + R*np.cos(theta)  # y-position of the center
#        print('xR:',xR)
#        print('________')
#        print('yR:',yR)
    if np.pi/2<azim<=np.pi:
        theta = (360*np.pi/180) - azim + np.pi/2  # rotation angle
#        print('theta angle;', theta)
        xR = X + R*np.sin(theta)  # x-position of the center
        yR = Y + R*np.cos(theta)  # y-position of the center
#        print('xR:',xR)
#        print('________')
#        print('yR:',yR)
    if np.pi<azim<=3*np.pi/2:
        theta = 2*np.pi - azim + np.pi/2  # rotation angle
#        print('theta angle;', theta)
        xR = X + R*np.sin(theta)  # x-position of the center
        yR = Y + R*np.cos(theta)  # y-position of the center
#        print('xR:',xR)
#        print('________')
#        print('yR:',yR)
    if 3*np.pi/2<azim<=2*np.pi: 
        theta = 2*np.pi - azim + np.pi/2  # rotation angle
#        print('theta angle;', theta)
        xR = X + R*np.sin(theta)  # x-position of the center
        yR = Y + R*np.cos(theta)  # y-position of the center
#        print('xR:',xR)
#        print('________')
#        print('yR:',yR)

    t = np.linspace(0, 2*np.pi, 100)

    x = xR + a*np.cos(azim)*np.cos(t) - b*np.sin(azim)*np.sin(t)
    y = yR + a*np.sin(azim)*np.cos(t) + b*np.cos(azim)*np.sin(t)
  #  plt.axis('equal')
  
    for i in range(len(x)):
        z = 0
        x[i],y[i], z = CarttoGeo(x[i],y[i],Z)
        print("x:",x[i])
        print('y:',y[i])
        
    plt.plot(x, y, color=color)  # rotated ellipse
    # plt.fill_between(x,y, color=color, alpha=0.4)

    return

###############################################################################

def Center(a, b, R, color):
    azim = np.linspace(0, 2*np.pi, 50)    
    for angle in azim:
        if 0<=angle<=np.pi/2:
            theta = (90*np.pi/180) - azim  # rotation angle
    #        print('theta angle;', theta)
            xR = R*np.sin(theta)  # x-position of the center
            yR = R*np.cos(theta)  # y-position of the center
    #        print('xR:',xR)
    #        print('________')
    #        print('yR:',yR)
        if np.pi/2<angle<=np.pi:
            theta = (360*np.pi/180) - azim + np.pi/2  # rotation angle
    #        print('theta angle;', theta)
            xR = R*np.sin(theta)  # x-position of the center
            yR = R*np.cos(theta)  # y-position of the center
    #        print('xR:',xR)
    #        print('________')
    #        print('yR:',yR)
        if np.pi<angle<=3*np.pi/2:
            theta = 2*np.pi - azim + np.pi/2  # rotation angle
    #        print('theta angle;', theta)
            xR = R*np.sin(theta)  # x-position of the center
            yR = R*np.cos(theta)  # y-position of the center
    #        print('xR:',xR)
    #        print('________')
    #        print('yR:',yR)
        if 3*np.pi/2<angle<=2*np.pi: 
            theta = 2*np.pi - azim + np.pi/2  # rotation angle
    #        print('theta angle;', theta)
            xR = R*np.sin(theta)  # x-position of the center
            yR = R*np.cos(theta)  # y-position of the center
    #        print('xR:',xR)
    #        print('________')
    #        print('yR:',yR)
            plt.scatter(xR,yR, color=color)
    return
            
# geo.overijssel.nl to look for db of rivers
# rijkswaterstaa.nl to look for bd of rivers

