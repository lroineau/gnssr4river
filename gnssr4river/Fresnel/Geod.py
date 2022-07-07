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

# rot = np.array([       -np.sin(lambda1),               np.cos(lambda1),              0     ],
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


#def sp3toGeo(sp3,x,y,z):
    """
    This function takes the coordinates given in a sp3 file and gives the elevation and azimut of the satellite seen from the receiver.
    
    Parameters
    ----------
    file: DataFrame
        sp3 file 
    x,y,z: cartesian coordinates of the receiver in meters
        sp3 file 
        
    Return
    ------
    azim, elev: float
        elevation an azimut of the satellite in degrees 
        
    """
 """   # Take the coordinates from the sp3, and set unit to meters
    xsat = sp3["x"] * 1000
    ysat = sp3["y"] * 1000
    zsat = sp3["z"] * 1000
    
    evel_sat = []
    azim_sat = []
    B = []
    
    for i in range(len(xsat)):
        X = xsat[i]-x
        Y = ysat[i]-y
        Z = zsat[i]-z
        
        B.append([X,Y,Z])
        B = np.array(B)
        
        rot = np.array([-np.sin(lambda1),  -np.sin(phi1)*np.cos(lambda1),   np.cos(phi1)*np.cos(lambda1)],
                       [ np.cos(lambda1),  -np.sin(phi1)*np.sin(lambda1),   np.cos(phi1)*np.sin(lambda1)],
                       [        0,                  np.cos(phi1),                   np.sin(phi1)        ])
        
        N = np.array([N1*np.cos(lambda1)*np.cos(phi1)],
                     [N1*np.sin(lambda1)*np.cos(phi1)],
                     [N1*(1-e**2)*np.sin(phi1)       ])
        
        A = np.dot(rot,B)
    
        C = np.add(N,A)
        
    return B
"""

def up(lon,lat):
    """
    This function gives the normal, east and north vectors for azimut calculs.
    
    Parameters
    ----------
    lat,lon: float
        location of latitude and longitudee in degrees 
        
    Return
    ------
    v_norm, East, North: array
        up unit vector, and local east and north unit vectors
    """
    # set unit to radians
    lon = m.radians(lon)
    lat = m.radians(lat)
    
    # up vector
    xo = np.cos(lon)*np.cos(lat)
    yo = np.sin(lon)*np.cos(lat)
    zo = np.sin(lat)
    v_norm = np.array([xo,yo,zo])    

    # north and east vectors
    North = np.zeros(3)
    East = np.zeros(3)
    North[0] = -np.sin(lat)*np.cos(lon)
    North[1] = -np.sin(lat)*np.sin(lon)
    North[2] = np.cos(lat)
    East[0] = -np.sin(lon)
    East[1] = np.cos(lon)
    East[2] = 0
    return v_norm, East, North


def elev_angle(v_norm, VSat_Rec):
    """
    This function gives the elevation of the satellite from the receiver.
    
    Parameters
    ----------
    v_norm: array
        unit vector in up direction
    Vsat_Rec: array
        Cartesian vector that points from receiver to the satellite in meters
        
    Return
    ------
    elev: float
        elevation angle in degrees
    """
    ang = np.arccos(np.dot(VSat_Rec,v_norm) / (norm(VSat_Rec)))
    angle = np.pi/2.0 - ang
    
    return angle*180/np.pi