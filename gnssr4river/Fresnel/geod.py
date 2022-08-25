# -*- coding: utf-8 -*-
"""
This code contains usefull functions for geodesic calculs.
@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import numpy as np
import math as m
import pandas as pd

# Some parameters of WGS84 projection...
a = 6378137.0             # Earth radius (m)
f = 1./298.257223563      # flattening factor
e = np.sqrt(2*f-f**2)     # first eccentricity
b = a*(1-f)
e_pr = np.sqrt((a**2 - b**2) / b**2) # second excentricity

###############################################################################

def mirrorMat(axis):
    """
    Constructs a transformation axis which reverses the coordinates of a given axis.
    
    Parameters
    ----------
    axis: int
        Either 0,1,2 for x,y,z respectively.
    
    Return
    ------
    mat: np.array
        A diagonal numpy matrix.
    """
    mat=np.identity(3)
    mat[axis,axis]=-1
    return mat

###############################################################################

def rotMat(axis,angle):
    """
    Constructs a transformation axis which rotates around a given axis by an angle
     
    Parameters
    ----------
    axis: int
        Either 0,1,2 for x,y,z respectively.
    angle: float 
        angle in radians.
    
    Return
    ------
    mat: np.array
        A diagonal numpy matrix.
    """
    cosang=np.cos(angle)
    sinang=np.sin(angle)
    mat=np.identity(3)
    idx=[i for i in range(3) if i != axis]
    mat[idx,idx]=[cosang,cosang]

    if axis%2 == 0:
        mat[idx[0],idx[1]]=sinang
        mat[idx[1],idx[0]]=-sinang
    else:
        mat[idx[0],idx[1]]=-sinang
        mat[idx[1],idx[0]]=sinang

    return mat

###############################################################################

def geo2cart(lon, lat, h):
    """
    Transform geographical coordinates (WGS84) to cartesian coordinates (WGS84).
    
    Parameters
    ----------
    lat,lon: float
        location of latitude and longitudee in degrees. 
    h: meters
        ellipsoïdal height.
        
    Return
    ------
    X,Y,Z: float
        cartesian coordinates.
    """
    # Put units to radians
    lon = np.radians(lon)
    lat = np.radians(lat)
    
    # Intermediate parameters
    slat = np.sin(lat)
    clat = np.cos(lat)
    N = a/(np.sqrt(1 - e**2*slat*slat))

    # Coordinates
    X = (N + h)*clat*np.cos(lon)
    Y = (N + h)*clat*np.sin(lon)
    Z = (N*(1 - e**2) + h)*slat
    
    return X,Y,Z

###############################################################################

def cart2geo(X,Y,Z):
    """
    This function transform cartesian coordinates to geographical coordinates (WGS84).
    
    Parameters
    ----------
    X,Y,Z: float
        cartesian coordinates.
        
    Return
    ------
    lat,lon: float
        location of latitude and longitudee in degrees.
    h: meters
        ellipsoïdal height.
    """
    lon = np.arctan2(Y, X)
    p = np.sqrt(X**2+Y**2)
    lat0 = np.arctan((Z/p)/(1-e**2))
    error = 1
    i=0 # make sure it doesn't go forever
    tol = 1e-10
    while error > tol and i < 6:
        n = a**2/np.sqrt(a**2*np.cos(lat0)**2+b**2*np.sin(lat0)**2)
        h = p/np.cos(lat0)-n
        lat = np.arctan((Z/p)/(1-e**2*n/(n+h)))
        error = np.abs(lat-lat0)
        lat0 = lat
        i+=1
    
    return np.degrees(lon), np.degrees(lat), h

###############################################################################

def elevazim(df_sp3,lon,lat,h):
    """
    This function takes the coordinates given in a sp3 file and gives the elevation 
    and azimut of the satellite seen from the receiver. Only visible satellites are
    returned.
    
    Parameters
    ----------
    df_sp3: DataFrame
        sp3 file.
    lon,lat: float 
        Longitude and latitude of receivers in degrees.
    h: float
        Heigth of the receiver in meters.
        
    Return
    ------
    df: DataFrame
        Elevation and azimuth of visible satellites in degrees.
    """
    # Take the coordinates from the sp3, and set unit to meters
    xsat = df_sp3["x"] * 1000
    ysat = df_sp3["y"] * 1000
    zsat = df_sp3["z"] * 1000
    # Add the other columns
    prn = df_sp3["prn"]
    week = df_sp3["week"]
    date = df_sp3["date"]
    tow = df_sp3["tow"]
    clock = df_sp3["clock"]
    sys = df_sp3["system"]

    # Get x,y,z of receiver
    x,y,z = geo2cart(lon,lat,h)
    # Make the output DataFrame
    a_list = list(range(1, len(xsat)))
    df = pd.DataFrame(columns = ['PRN','system','date','week', 'tow','clock', 'elevation', 'azimuth'], index = a_list)

    # Put unit to radians
    lon = m.radians(lon)
    lat = m.radians(lat)
    colat=np.pi/2-lat 
    #cosntruct transformation matrix to transform fro global to local NEH frame
    Rtrans=np.matmul(mirrorMat(0),np.matmul(rotMat(colat,1),rotMat(lon,2)))

    ## Calcul elevation and azimut for each satellite of the file ##

    # Vector between satellite and antenna
    X = xsat-x
    Y = ysat-y
    Z = zsat-z
        
    # Rotation matrix
    N,E,H=np.matmul(Rtrans,np.array([X,Y,Z]))

    # Calculation of the azimuth
    az = np.arctan2(E, N)*180/np.pi        
    # if it is negatif, we add 2pi
    for i in range(len(az)):
        if az[i] < 0:
            az[i] = 360 + az[i]
  
    # Calculation of the elevation
    rxy = np.sqrt(N**2 + E**2)
    el = (np.arctan2(H,rxy))*180/np.pi

    d = {'PRN':prn, 'week':week, 'date':date, 'tow':tow, 'clock':clock, 'system':sys, 'elevation':np.round(el),'azimuth':np.round(az)}
    df = pd.DataFrame(data=d)

    # Drop potential nan inside dataframe
    df = df.dropna()
    
    # Set type, overwise it is String
    df = df.astype({'week':'int','system':'str','tow':'float','clock':'float','azimuth':'int', 'elevation':'int','PRN':'int','date':'datetime64'})
    
    # Remove all values bellow 0° elevation, i.e satellites that are not visible
    df = df[df['elevation'] > 0]
    
    return df