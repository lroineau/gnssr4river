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

# Project import
from FresnelZone import FirstFresnelZone, plotEllipse 
from IterFresnel import *

# Usefull constants
c = 299792458             # m.s-1 Speed of light
L1_freq = 1575.42e6       # Hz L1 frequency
L2_freq = 1227.60e6       # Hz L2 frequency
lambda_L1 = c/L1_freq     # m wavelenght for L1
lambda_L2 = c/L2_freq     # m wavelenght for L2

###############################################################################

def image_spoof(self, tile):
    '''
    This function reformats web requests from OSM for cartopy
    Heavily based on code by Joshua Hrisko at:
    https://makersportal.com/blog/2020/4/24/geographic-visualizations-in-python-with-cartopy
    '''

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
        location of latitude and longitudee in degrees        

    Return
    ------
    dist: float
        distance to edge from centre in meters

    """
    dist_cnr = np.sqrt(2*dist**2)
    top_left = cgeo.Geodesic().direct(points=(lon,lat),azimuths=-45,distances=dist_cnr)[:,0:2][0]
    bot_right = cgeo.Geodesic().direct(points=(lon,lat),azimuths=135,distances=dist_cnr)[:,0:2][0]

    extent = [top_left[0], bot_right[0], bot_right[1], top_left[1]]

    return extent

###############################################################################

def osm_image(lon,lat,dist,sitename=None,style='satellite',save=False):
    """
    This function makes OpenStreetMap satellite or map image.    
 
    Parameters
    ----------
    lat,lon: float
        location of latitude and longitudee in degrees   
    dist: float
        distance to edge from centre in meters  
    sitename: String
        name of the location to be written on plot
    style: String
        decide wether background is "map" or "satellite"
    save: boolean
        decide wether the return image will be saved or not
     
    Return
    ------
    Plot with map or satellite background
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

    ax.set_title(f'{sitename} ({lat},{lon})',fontsize=15)

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
            print("Directory {} created".format(dirName))
        else:    
            print("Directory {} already exists".format(dirName))

        fig.savefig(f'Img/{sitename}_{style}_r{dist}_scale{scale}.jpg', dpi=150, bbox_inches='tight')

    return


###############################################################################
# For test with Dinkel coordinates
# osm_image(6.951909549393194,52.42752035169519,250,sitename='TestSat',style='satellite',save=False)
# osm_image(6.951909549393194,52.42752035169519,250,sitename='TestMap',style='map',save=False)