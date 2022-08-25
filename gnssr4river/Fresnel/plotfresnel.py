# -*- coding: utf-8 -*-
"""
This code is to plot the reflection points on a map for vizualisation purpose

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Basics import
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.geodesic as cgeo
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import io
from urllib.request import urlopen, Request
from PIL import Image
from pathlib import Path

# Project import
from fresnelzone import *
from iterfresnel import *
from intersect import *

###############################################################################

def image_spoof(self, tile):
    """
    This function reformats web requests from OSM for cartopy.
    
    Heavily based on code by Joshua Hrisko at:
    https://makersportal.com/blog/2020/4/24/geographic-visualizations-in-python-with-cartopy
    """    
    url = self._image_url(tile)                # get the url of the street map API
    req = Request(url)                         # start request
    req.add_header('User-agent','Anaconda 3')  # add user agent to request
    fh = urlopen(req) 
    im_data = io.BytesIO(fh.read())            # get image
    fh.close()                                 # close url
    img = Image.open(im_data)                  # open image with PIL
    img = img.convert(self.desired_tile_form)  # set image format
    
    return img, self.tileextent(tile), 'lower' # reformat for cartopy

###############################################################################

def calc_extent(lon,lat,dist):
    """
    This function calculate extent of map
    
    Parameters
    ----------
    lat,lon: float
        Location of latitude and longitudee in degrees.   

    Return
    ------
    dist: float
        Distance to edge from centre in meters.
    """
    dist_cnr = np.sqrt(2*dist**2)
    top_left = cgeo.Geodesic().direct(points=(lon,lat),azimuths=-45,distances=dist_cnr)[:,0:2][0]
    bot_right = cgeo.Geodesic().direct(points=(lon,lat),azimuths=135,distances=dist_cnr)[:,0:2][0]

    extent = [top_left[0], bot_right[0], bot_right[1], top_left[1]]

    return extent

###############################################################################

def osm_image(lon,lat,dist,style='satellite',sitename=None,shp=True,save=False):
    """
    This function makes OpenStreetMap satellite or map image.    
 
    Parameters
    ----------
    lat,lon: float
        Location of latitude and longitudee in degrees.
    dist: float
        Distance to edge from centre in meters.
    style: String
        Decide wether background is "map" or "satellite".
    sitename: String
        Name of the location to be written on plot.
    save: boolean
        Decide wether the return image will be saved or not.
     
    Return
    ------
    Plot with map or satellite background.
    """
    # Check which background to apply
    if style=='map':
        ## MAP STYLE
        cimgt.OSM.get_image = image_spoof # reformat web request for street map spoofing
        img = cimgt.OSM() # spoofed, downloaded street map
    elif style =='satellite':
        # SATELLITE STYLE
        cimgt.QuadtreeTiles.get_image = image_spoof # reformat web request for street map spoofing
        img = cimgt.QuadtreeTiles() # spoofed, downloaded street map
    else:
        print('no valid style')

    plt.close('all')
    fig = plt.figure(figsize=(10,10)) # open matplotlib figure
    ax = plt.axes(projection=img.crs) # project using coordinate reference system (CRS) of street map
    data_crs = ccrs.PlateCarree()
    
    ax.set_title(f'{style} overview of {sitename} \n ({lat},{lon})',fontsize=15)

    # auto-calculate scale
    scale = int(120/np.log(dist))
    scale = (scale<20) and scale or 19

    extent = calc_extent(lon,lat,dist*2.5)
    ax.set_extent(extent) # set extents
    ax.add_image(img, int(scale)) # add OSM with zoom specification

    gl = ax.gridlines(draw_labels=True, crs=data_crs,
                        color='k',lw=0.5)

    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    
    if save==True:
        # Create directory
        dirName = 'Img'
        # Create target Directory if don't exist
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print(f"Image stored in {dirName}")
        else:    
            print(f"Image stored in {dirName}")

        fig.savefig(f'Img/{sitename}_{style}_scale{scale}.jpg', dpi=150, bbox_inches='tight')

    return

###############################################################################

def fresnelMap(lon,lat,dist=25,filezone=None,filewb=None,azmin=None,azmax=None,alti=False,save=False,sitename=None,dirName=None):
    """
    This function add Fresnel zones to a satellite background.
 
    Parameters
    ----------
    lat,lon: float
        Location of latitude and longitudee in degrees.
    dist: float
        Distance to edge from centre in meters, value 25 is default.
    filezone: str
        Path of the file for the Fresnel Zone, can be either gpkg or shp or kml.
    filewb: str
        Path of the file for the waterbody, can be either gpkg or shp or kml.
    azmin,azmax: int
        Define custom azimuth mask, min can be superior to max if mask overlap azimuth 0.
    alti: Boolean
        If True, return information on altimeter tracks.
    save: Boolean
        Decide wether the return image will be saved or not.
    sitename: str
        Name of the location, used for plot title and name of file if saved.
    dirName: str
        Directory where files will be stored.
        
    Return
    ------
    Plot of Fresnel Zones with satellite background.
    """
    cimgt.QuadtreeTiles.get_image = image_spoof # reformat web request for street map spoofing
    img = cimgt.QuadtreeTiles() # spoofed, downloaded street map

    plt.close('all')
    fig = plt.figure(figsize=(10,10)) # open matplotlib figure
    ax = plt.axes(projection=img.crs) # project using coordinate reference system (CRS) of street map
    data_crs = ccrs.PlateCarree()
    
    # Handling Fresnel Zones
    if filezone==None:
        raise Exception("No geo file to plot given")  
    else:
        Fzone = gpd.read_file(filezone)
        FzoneM=Fzone.to_crs(crs=3857)
        filename = Path(filezone).stem
        
    # Handling Water Body
    if filewb==None:
        print("No file for waterbody given")
    else:
        wb = gpd.read_file(filewb)
        ax.add_geometries(wb.geometry,crs = ccrs.PlateCarree(),facecolor='dodgerblue', edgecolor='dodgerblue',alpha=0.75)
        Fzone, s, waz=intersect_wb(Fzone, wb)
    
    # Handling Altimeter tracks
    if alti==True:
        if filewb==None:
            intersect_altitrack(FzoneM)
        else:
            intersect_altitrack(FzoneM,filewb=filewb)

    # If user wants its own azimuth mask
    if azmin!=None and azmax!=None:
        Fzone = Fzone.loc[((Fzone.azimuth > azmin) & (Fzone.azimuth < azmax))]
            
    # List of elevation present in DataFrame
    d = []
    for i in range(len(Fzone.elevation.unique())):
        d.append(Fzone[Fzone.elevation == Fzone.elevation.unique()[i]])

    # List of color
    c=['yellow','blue','red','green']
    z=5
    # Plot Fresnel Zones by elevation
    for i in range(len(d)):
        ax.add_geometries(d[i].geometry,crs = ccrs.PlateCarree(),facecolor=c[i],edgecolor=c[i],alpha=0.35,zorder=z)
        z+=5
    n_elev = Fzone.elevation.unique()
    n_elev = n_elev.tolist()
    co = []
    for i in range(len(n_elev)):
        co.append(c[i])
    ax.set_title(f'Fresnel zones for {n_elev}° elevation, respectively {co}',fontsize=20,pad=10)  
    
    # auto-calculate scale
    scale = int(120/np.log(dist))
    scale = (scale<20) and scale or 19

    extent = calc_extent(lon,lat,dist*2.5)
    ax.set_extent(extent) # set extents
    ax.add_image(img, int(scale)) # add OSM with zoom specification

    gl = ax.gridlines(draw_labels=True, crs=data_crs,
                        color='k',lw=0.5)

    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    
    # Handling save
    if save==True:
        # Create directory
        if dirName==None:
            dirName = 'Img'
        # Create target Directory if doesn't exist
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Image will be stored in {dirName}")
        if sitename!=None:
            fig.savefig(f'Img/{sitename}_{filename}.jpg', dpi=150, bbox_inches='tight')
        else:
            fig.savefig(f'Img/{filename}.jpg', dpi=150, bbox_inches='tight')            

    if 's' and 'waz' in vars():
        return s,waz
    else:
        return