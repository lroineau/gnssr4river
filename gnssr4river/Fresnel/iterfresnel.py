# -*- coding: utf-8 -*-
"""
To iterate calculation of Fresnel Zone and show reflexion points for GNSS-R

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede, NL), Aug 26, 2022
"""


# Usefull librairies
import numpy as np
import matplotlib.pyplot as plt
from private.FresnelZone import FirstFresnelZone, plotEllipse
from private.GetOrbits import *
from private.PlotFresnel import *
import geopandas as gpd

# Usefull constants
c = 299792458             # m.s-1 Speed of light
L1_GPS = 1575.42e6         # Hz L1 frequency for GPS
L2_GPS = 1227.60e6         # Hz L2 frequency for GPS
L1_Glo = 1575.42e6         # Hz L1 frequency for Glonass
L2_Glo = 1227.60e6         # Hz L2 frequency for Glonass
lambda_L1_GPS = (c/L1_GPS)     # m wavelenght for L1 GPS
lambda_L2_GPS = (c/L2_GPS)     # m wavelenght for L2 GPS
lambda_L1_Glo = (c/L1_Glo)     # m wavelenght for L1 Glonass
lambda_L2_Glo = (c/L2_Glo)     # m wavelenght for L2 Glonass
crs = {'init': 'epsg:4326'}
###############################################################################

def IterFresnel(freq, h, elev_list, lon, lat, azim_min, azim_max, df_sp3, nb=360):
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
    lon,lat: float
        position of the receiver in geographical coordinates (degrees)
    h: float
        hight of the receiver in meters
    azim_min, azim_max: float
        minimum and maximum azimut for calcule  
    nb: int
        number of ellipse to plot
    df_sp3: DataFrame
        unsorted sp3 file of satellites.
    
    Return
    ------
    A plot with all Fresnel zones
    """
    # Set azimut and color
    # azim = np.linspace(azim_min, azim_max, nb)
    color = ['yellow','blue','red','green']

    # Empty list to store Polygon and Polygons area
    pol = []    
    # l_area = []
    
    # Check if elevation is a list or a single value
    if isinstance(elev_list, list):
        print('More than one elevation given, processing...')
        # Calcule First Fresnel Surface for each elevation angles
        c = 0
        for elev in elev_list:
            df_int = df_sp3.loc[(df_sp3['elevation'] == elev)]
            azim = df_int["azimuth"]
            azim = azim.sort_values()
            # Prevent to calculate several time the same azimuth
            list_azim = []
            for val in azim:
                if val not in list_azim:
                    list_azim.append(val)
            for angle in list_azim:
                a,b,R = FirstFresnelZone(freq, h, elev)
                p = plotEllipse(a,b,R,lon,lat,h,angle,color[c])
                pol.append(p)
                # l_area.append(area)
            c+=1
                    
    else:
        print('only one elevation given, processing...')
        elev = elev_list
        df_sp3 = df_sp3.loc[(df_sp3['elevation'] == elev)]
        azim = df_sp3["azimuth"]
        azim = azim.sort_values()
        # Prevent to calculate several time the same azimuth
        list_azim = []
        for val in azim:
            if val not in list_azim:
                list_azim.append(val)
        # Calcule First Fresnel Surface for each elevation angle
        for angle in list_azim:
            a,b,R = FirstFresnelZone(freq, h, elev)
            p = plotEllipse(a,b,R,lon,lat,h,angle, 'green')
            pol.append(p)
            # l_area.append(area)

    # Shp with GeoPandas
    gdf = gpd.GeoDataFrame(crs=crs,geometry=pol)
    #gdf.assign(area=l_area)
    gdf.to_file(f"Fresnel{elev_list}_{h}.shp", driver='ESRI Shapefile')

    return gdf       


# For test with Dinkel coordinates
# IterFresnel(L1_freq, 2, 25, 0, 360, nb=50)
# IterFresnel(L1_freq, 2, [5,10,15], 0, 360, nb=50)