# -*- coding: utf-8 -*-
"""
This code contains usefull functions for geodesic calculs.

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt

# Some parameters of WGS84 projection...
a = 6378137.0             # Earth radius (m)
f = 1./298.257223563      # flattening factor
e = np.sqrt(2*f-f**2)     # first eccentricity
b = a*(1-f)
e_pr = np.sqrt((a**2 - b**2) / b**2) # second excentricity

###############################################################################

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



def elevazim(df_sp3,lon,lat,h):
    """
    This function takes the coordinates given in a sp3 file and gives the elevation and azimut of the satellite seen from the receiver.
    
    Parameters
    ----------
    df_sp3: DataFrame
        sp3 file 
    lon,lat: float 
        longitude and latitude of receivers in degrees
    h: float
        heigth of the receiver in meters
        
    Return
    ------
    df: DataFrame
        elevation an azimut of the satellite in degrees 
        
    """
    # Take the coordinates from the sp3, and set unit to meters
    xsat = df_sp3["x"] * 1000
    ysat = df_sp3["y"] * 1000
    zsat = df_sp3["z"] * 1000
    # Add the other columns
    prn = df_sp3["prn"]
    week = df_sp3["week"]
    tow = df_sp3["tow"]
    clock = df_sp3["clock"]

    print("Calculation of satellites coordinates in progress, may take a while...")

    # Get x,y,z of receiver
    x,y,z = GeotoCart(lon,lat,h)
    
    # Make the output DataFrame
    a_list = list(range(1, len(xsat)))
    df = pd.DataFrame(columns = ['PRN', 'week', 'tow','clock', 'elevation', 'azimuth'], index = a_list)

    # Put unit to radians
    lon = m.radians(lon)
    lat = m.radians(lat)

    ## Calcul elevation and azimut for each satellite of the file ##
    for i in range(len(xsat)):
        # Vector between satellite and antenna
        X = xsat[i]-x
        Y = ysat[i]-y
        Z = zsat[i]-z
        
        # Rotation matrix
        N = -np.sin(lat)*np.cos(lon)*X - np.sin(lat)*np.sin(lon)*Y + np.cos(lat)*Z
        E = - np.sin(lon)*X + np.cos(lon)*Y
        Z = np.cos(lat)*np.cos(lon)*X + np.cos(lat)*np.sin(lon)*Y + np.sin(lat)*Z
 
        # Calculation of the azimuth
        az = np.arctan2(E, N)*180/np.pi        
        # if it is negatif, we add 2pi
        if az < 0:
            az = 360 + az
  
        # Calculation of the elevation
        r = np.sqrt(X**2 + Y**2 + Z**2)
        # el = np.arcsin(Z/r)
        angle = np.arccos(Z/r)
        el2 = (np.pi/2.0 - angle)*180/np.pi

        # Putting everything in the DataFrame        
        sat = prn[i]
        we = week[i]
        tw = tow[i]
        cl = clock[i]
        df.loc[f'{i}'] = [f'{sat}',f'{we}',f'{tw}',f'{cl}',round(el2),round(az)] 
        
    # Drop potential nan inside dataframe
    df = df.dropna()
    
    # Set type to float, overwise it is String
    df = df.astype({'week':'float','tow':'float','clock':'float','azimuth':'float'})
    return df