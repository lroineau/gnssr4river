# -*- coding: utf-8 -*-
"""
This code is to retrieve the orbits of GPS and GLONASS satellites

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

import wget
import os
import numpy as np
import zlib
import pandas as pd
import datetime as d


# Some parameters of WGS84 projection...
a = 6378137.0               # Earth radius (m)
f = 1./298.257223563      # flattening factor
e = np.sqrt(2*f-f**2)       # first eccentricity


def gpsweek(year, month, day, hour, minute, second):
    """
    Return the gps week and secondes for a given date (adapted from a code by Kristin Larson). 

    Parameters
    ----------
    year, month, day, hour, minute, second : int
        Year (xxx), month, day, hour, minute, second to be change to gps week (default is current daytime)

    Return
    ------
    GPS_wk, GPS_sec_wk : int       
        the gps week and second of the week
    """
    UT=hour+minute/60.0 + second/3600. 
    if month > 2:
        y=year
        m=month
    else:
        y=year-1
        m=month+12
        
    JD=np.floor(365.25*y) + np.floor(30.6001*(m+1)) + day + (UT/24.0) + 1720981.5
    GPS_wk=np.floor((JD-2444244.5)/7.0);
    GPS_wk = np.int(GPS_wk)
    GPS_sec_wk=np.rint((((JD-2444244.5)/7)-GPS_wk)*7*24*3600)            
     
    return GPS_wk, GPS_sec_wk


def retrieve_orbits(year=None, month=None, day=None, hour=None, minute=None, second=None):
    """
    Retrieve the orbits of gnss constellations from ESA site. 

    Parameters
    ----------
    
    Return
    ------
    
    """
    # if the date is not given by the user, take current date
    if year==None:
        cur_date = d.datetime.now()
        cur_date = cur_date.strftime("%Y,%m,%d,%H,%M,%S")
        
        cur_date = cur_date.split(",")
        year, month, day, hour, minute, second = int(cur_date[0]), int(cur_date[1]), int(cur_date[2]), int(cur_date[3]), int(cur_date[4]), int(cur_date[5])
    
    # Retriecve first the gps week
    GPS_wk, GPS_sec_wk = gpsweek(year, month, day, hour, minute, second)
    
    # Create directory
    dirName = 'Orbits'
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory {} created".format(dirName))
    else:    
        print("Directory {} already exists".format(dirName))
        
    # data link
    url = 'http://navigation-office.esa.int/products/gnss-products/{}/esu22156_06.sp3.Z'.format(GPS_wk)
    # Check if file already exists 
    if os.path.exists('Orbits/esu22156_06.sp3.Z') is True:
        print("File already exists")
    # Else it is downloaded
    else:
        wget.download(url, out=dirName)
        if os.path.exists('Orbits/esu22156_06.sp3.Z') is True:
            print("Data download success")
        else:
            print("Fail to retrieve data")

    # uncompress file
    file = open('Orbits/esu22156_06.sp3.Z', 'rb').read()
    orb = zlib.decompress(file)
    
    return orb

     

def read_sp3(file):
    """
    Read the sp3 file and turn it to a DataFrame.

    Parameters
    ----------
    file : String
        the sp3 file
    
    Return
    ------
    df_sp3: DataFrame
        the content of the sp3 on a DataFrame
    """
    try:      
        f = open('Orbits/{}'.format(file))
        raw = f.read()
        f.close()
        lines  = raw.splitlines()
        nprn = int(lines[2].split()[1])
        lines  = raw.splitlines()[22:-1]
        epochs = lines[::(nprn+1)]
        nepoch =  len(lines[::(nprn+1)])
        week, tow, x, y, z, clock, prn = np.zeros((nepoch*nprn, 7)).T
        for i in range(nepoch):
            year, month, day, hour, minute, second = np.array(epochs[i].split()[1:], dtype=float)
            week[i*nprn:(i+1)*nprn], tow[i*nprn:(i+1)*nprn] = \
				gpsweek(year, month, day, hour, minute, second)
            for j in range(nprn):
                prn[i*nprn+j] =  int(lines[i*(nprn+1)+j+1][2:4])
                x[i*nprn+j] = float(lines[i*(nprn+1)+j+1][4:18])
                y[i*nprn+j] = float(lines[i*(nprn+1)+j+1][18:32])
                z[i*nprn+j] = float(lines[i*(nprn+1)+j+1][32:46])
                clock[i*nprn+j] = float(lines[(i)*(nprn+1)+j+1][46:60])
    except:
        print('sorry - the sp3file does not exist')
        week,tow,x,y,z,prn,clock=[0,0,0,0,0,0,0]
		
    df_sp3 = pd.DataFrame({"week":week,
                           "tow":tow, 
                           "x":x,
                           "y":y,
                           "z":z,
                           "prn":prn,
                           "clock":clock})
        
    return df_sp3