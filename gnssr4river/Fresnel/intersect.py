# -*- coding: utf-8 -*-
"""
Contains function for intersecting data to Fresnel zones

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Basics import

import geopandas as gpd
import pandas as pd
import numpy as np

def intersect_wb(Fzone, wb):
    """
    Check intersection between Fresnel and water body.

    Parameters
    ----------
    Fzone : GeoDataFrame
        GeoDataFrame of the Fresnel zones.
    wb : GeoDataFrame
        GeoDataFrame of the water body.

    Returns
    -------
    Fzone : GeoDataFrame
        GeoDataFrame of the Fresnel zones with only ellipses that intersect the
        water body.
    s : int
        Number of ellipses completely within the water body.
    waz : list
        Windows of azimuth for observations.
        
    """
    intersection_gdf = Fzone.overlay(wb, how='intersection')
    Fzone = Fzone[Fzone.ID.isin(intersection_gdf['ID'])]
    Az_zone = Fzone['azimuth']
    waz = window_az(Az_zone)
    FzoneWithin = gpd.sjoin(Fzone, wb, how='inner', op='within')
    s = len(FzoneWithin)
    print(f'Number of Fresnel zones completely within water body is {s}')
    return Fzone, s, waz


def intersect_altitrack(FzoneM,filewb=None):   
    """
    Check intersection between altimeter tracks from Sentinel 3 satellites.

    Parameters
    ----------
    Fzone : GeoDataFrame
        GeoDataFrame of the Fresnel zones, geometry in meters
    filewb : GeoDataFrame, optional
        GeoDataFrame of the water body. The default is None.

    Returns
    -------
    None.

    """
    # Sentinel 3
    gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
    alti3A = gpd.read_file('sentinel3/S3A_rel_orbit_ground_track_10sec_v1_4.kml', driver='KML')
    alti3B = gpd.read_file('sentinel3/orbits/S3B_rel_orbit_ground_track_10sec_v1_4.kml', driver='KML')
    alti3AM=alti3A.to_crs(crs=3857)
    alti3BM=alti3B.to_crs(crs=3857)
    if 'ax' in vars():
        ax.add_geometries(alti3A.geometry, crs = ccrs.PlateCarree(), facecolor='red', edgecolor='red')
        ax.add_geometries(alti3B.geometry, crs = ccrs.PlateCarree(), facecolor='red', edgecolor='green')
        
    # Check if alti tracks intersect water body
    if filewb==None:
        pass
    else:
        dwb = gpd.read_file(filewb)
        intersection_3A = alti3A.overlay(dwb, how='intersection')
        if len(intersection_3A.index) == 0:
            print("No intersection with Sentinel 3A tracks")
        else:
            print("Sentinel 3A tracks intersect the water body")
        intersection_3B = alti3B.overlay(dwb, how='intersection')
        if len(intersection_3B.index) == 0:
            print("No intersection with Sentinel 3B tracks")
        else:
            print("Sentinel 3B tracks intersect the water body")
          
    # Return minimal distance between the reflection points and the altimeter tracks
    min_dist3A = round(min(FzoneM.distance(alti3AM)))
    min_dist3B = round(min(FzoneM.distance(alti3BM)))
    print(f'Nearest altimeter track from Sentinel 3 satellites is at {min(min_dist3A,min_dist3B)/1000} km')
    print('________________________________\n')
    return

def window_az(az):
    """
    Returns list of azimuth window of observations, i.e. azimuth range where 
    there are more than 3 reflections. If two azimuths that follow one another
    have a difference greater than 15° degrees, they are considered to be part
    of 2 different windows.

    Parameters
    ----------
    az : DataFrame column
        DataFrame column of azimuths.

    Returns
    -------
    waz : list
        The different windows of observation.
        
    """
    az = az.unique()
    #az = az.to_numpy()
    az.sort()
    
    # First check number of windows
    m = [[az[0]]]
    res=[[az[0]]]
    for x in range(len(az[1:])):
        if az[x+1] - m[-1][0] < 15:
            res[-1].append(az[x+1])
        else:
            res.append([az[x+1]])
        m=[[az[x+1]]]

    # Remove windows with less than 3 values
    a=[]
    for i in range(len(res)):
        if len(res[i])>=3:
            a.append(res[i])
    
    # Get min and max for each window
    azmin=[]
    azmax=[]
    for i in range(len(a)):
        azmin.append(min(a[i]))
        azmax.append(max(a[i]))
    
    rev=False    
    # Check if last window doesn't overlap first one (because of circle)
    for i in range(len(azmin)):
        for j in range(len(azmax)):
            if azmax[j]>350 and azmin[i]<10:
                azmax[j]=azmax[0]
                rev=True
    
    # Get final windows
    waz=[]
    print('Window of observation:')   
    if rev==True:
        for d in range(1,len(azmin)):
            print('Azimuth between '+str(azmin[d])+'° and '+str(azmax[d])+'°')
            waz.append((azmin[d],azmax[d]))
    else:
        for d in range(len(azmin)):
            print('Azimuth between '+str(azmin[d])+'° and '+str(azmax[d])+'°')
            waz.append((azmin[d],azmax[d]))

    return waz