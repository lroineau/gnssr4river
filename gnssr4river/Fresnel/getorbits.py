# -*- coding: utf-8 -*-
"""
This code is to retrieve the orbits of GPS and GLONASS satellites
@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import requests
import os
import numpy as np
import unlzw3
from pathlib import Path
import pandas as pd
from datetime import datetime 
from io import BytesIO
from geod import *
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

###############################################################################

def gpsweek(date=datetime.now()):
    """
    Return the gps week and secondes for a given date (adapted from a code by Kristin Larson). 
    Parameters
    ----------
    date: datetime
        date to use for determining the gps week
    Return
    ------
    GPS_wk, GPS_sec_wk : int       
        the gps week and second of the week
    """

    hour=date.hour
    minute=date.minute
    second=date.second
    month=date.month
    year=date.year
    day=date.day
    UT=hour+minute/60.0 + second/3600. 
    if month > 2:
        y=year
        m=month
    else:
        y=year-1
        m=month+12
        
    JD=np.floor(365.25*y) + np.floor(30.6001*(m+1)) + day + (UT/24.0) + 1720981.5
    GPS_wk=np.floor((JD-2444244.5)/7.0);
    GPS_wk = int(GPS_wk)
    GPS_sec_wk=np.rint((((JD-2444244.5)/7)-GPS_wk)*7*24*3600)            
     
    return GPS_wk, GPS_sec_wk

###############################################################################

def dateUTC(GPS_wk, GPS_sec_wk):
    """
    Inverse function that gives the UTC date from data contain in a GNSS file    
    Parameters
    ----------
    date: datetime
        date to use for determining the gps week
    Return
    ------
    GPS_wk, GPS_sec_wk : int       
        the gps week and second of the week
    """


    hour=date.hour
    minute=date.minute
    second=date.second
    month=date.month
    year=date.year
    day=date.day
    UT=hour+minute/60.0 + second/3600. 
    if month > 2:
        y=year
        m=month
    else:
        y=year-1
        m=month+12
        
    JD=np.floor(365.25*y) + np.floor(30.6001*(m+1)) + day + (UT/24.0) + 1720981.5
    GPS_wk=np.floor((JD-2444244.5)/7.0);
    GPS_wk = int(GPS_wk)
    GPS_sec_wk=np.rint((((JD-2444244.5)/7)-GPS_wk)*7*24*3600)            
     
    return GPS_wk, GPS_sec_wk

###############################################################################

def retrieve_orbits(date=datetime.now()):
    """
    Retrieve the orbits of GPS and Glonass constellations from ESA site. If date is not specified, takes current date.
    Parameters
    ----------
    date: datetime
        The date used to determine the GPS week for which the orbits will be downloaded (default takes the current date)  
    Return
    ------
    
    """
    # if user gives an incorrect date
    if date > datetime.now():
        raise RuntimeError("Cannot give a date that is yet to come !")  

    # Retrieve first the gps week
    GPS_wk, GPS_sec_wk = gpsweek(date)    
    print('GPS week is:', GPS_wk)
    
    # Create directory
    dirName = 'Orbits'
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory {} created".format(dirName))
    else:    
        print("Directory {} already exists".format(dirName))
        
    # Name of file to retrieve
    t = str(int(GPS_sec_wk/86400))
    # Zipped name
    filenameZ = f'esu{GPS_wk}' + t + '_00.sp3.Z' 
    print('filename is:',filenameZ)
    # Unzipped name    
    filename = f'esu{GPS_wk}' + t + '_00.sp3'         
    # data link   
    url = f'http://navigation-office.esa.int/products/gnss-products/{GPS_wk}/{filenameZ}'
    # Check if file already exists 
    if os.path.exists('Orbits/{}'.format(filenameZ)) is True:
        print("File already exists")
    # Else it is downloaded
    else:
        r = requests.get(url)
        open(dirName +'/' +filenameZ, 'wb').write(r.content)
        if os.path.exists('Orbits/{}'.format(filenameZ)) is True:
            print("Data download success")
        else:
            print("Fail to retrieve data")
        
    orb = unlzw3.unlzw(Path(f'Orbits/esu{GPS_wk}{t}_00.sp3.Z'))

    return orb

###############################################################################  

def read_sp3(file):
    """
    Read the sp3 file and turn it to a DataFrame. Can also be a Bytes file if unzipped directly with unlzw3.
    
    Parameters
    ----------
    file : String or Bytes
        the sp3 file or the local unzipped file in Python
    
    Return
    ------
    df_sp3: DataFrame
        the content of the sp3 on a DataFrame
    """
    # See if the file is a sp3 on hard drive or directly a Bytes in Python
    if type(file)==bytes:        
        f = BytesIO(file) 
    else:
        f = open('Orbits/{}'.format(file))
    
    # Read the file 
    try:      
        raw = f.read()
        f.close()
        lines  = raw.splitlines()
        nprn = int(lines[2].split()[1])
        lines  = raw.splitlines()[22:-1]
        epochs = lines[::(nprn+1)]
        nepoch =  len(lines[::(nprn+1)])
        week, tow, x, y, z, clock, prn = np.zeros((nepoch*nprn, 7)).T
        sys = []
        for i in range(nepoch):
            year, month, day, hour, minute, second = np.array(epochs[i].split()[1:], dtype=float)
            dtepoch=datetime(year=int(year),month=int(month),day=int(day),hour=int(hour),second=int(second))
            week[i*nprn:(i+1)*nprn], tow[i*nprn:(i+1)*nprn] = \
				gpsweek(dtepoch)
            for j in range(nprn):
                if str(lines[i*(nprn+1)+j+1][1:2])=="b'R'":
                    sys.append('GLONASS')
                elif str(lines[i*(nprn+1)+j+1][1:2])=="b'G'":
                    sys.append('GPS')
                prn[i*nprn+j] =  int(lines[i*(nprn+1)+j+1][2:4])
                x[i*nprn+j] = float(lines[i*(nprn+1)+j+1][4:18])
                y[i*nprn+j] = float(lines[i*(nprn+1)+j+1][18:32])
                z[i*nprn+j] = float(lines[i*(nprn+1)+j+1][32:46])
                clock[i*nprn+j] = float(lines[(i)*(nprn+1)+j+1][46:60])
                
    # if file not found
    except Exception as exc:
        print('sorry - the sp3file does not exist')
        week,tow,x,y,z,prn,clock=[0,0,0,0,0,0,0]
		
    # Set the DataFrame
    df_sp3 = pd.DataFrame({"date":dtepoch,
                           "system":sys,
                           "week":week.astype(int),
                           "tow":tow, 
                           "x":x,
                           "y":y,
                           "z":z,
                           "prn":prn.astype(int),
                           "clock":clock})
        
    return df_sp3

###############################################################################

def elevationSort(df_sp3, elev):
    """
    This function takes the sp3 DataFrame of satellites, and keep only good values of elevation.
    
    Parameters
    ----------
    df_sp3: DataFrame
        sp3 DataFrame of the satellite.
    elev: float or list
        elevation in degrees, should be one number or a list of 2 to give range.
    Returns
    -------
    df_sort: DataFrame
        DataFrame with wanted value for Fresnel.
    """
    # First remove all elevation bellow 0
    df_sp3 = df_sp3[df_sp3['elevation'] > 0]
    
    # If elev is a list
    if isinstance(elev, list):
        l = len(elev)
        print("Input is a list")
        
        elev_min = elev[0] - 1
        elev_max = elev[l-1] + 1
        
        print("Minimum elevation is:", elev_min, "\nMaximum elevation is:", elev_max)
        
        # Just check if user didn't miss input
        if elev_min>elev_max:
            raise Exception("Minimum elevation greater than maximum elevation")  

        df_sp3 = df_sp3.loc[(df_sp3['elevation'] >= elev_min) & (df_sp3['elevation'] <= elev_max)]

    # If elev is only one number
    else:
        
        print("Input is a number")
        
        # Take some range to not remove all elevation values
        elev_min = elev - 1
        elev_max = elev + 1
        
        print("Minimum elevation is:", elev_min, "\nMaximum elevation is:", elev_max)
        
        df_sp3 = df_sp3.loc[(df_sp3['elevation'] >= elev_min) & (df_sp3['elevation'] <= elev_max)]
    
    return df_sp3

###############################################################################

def clean_dir():
    """
    Simple function to clean directory where sp3 files arre stored.
    """
    dir = 'Orbits'
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))
    return

###############################################################################

def skyplot(df,GPS=True,GLONASS=True):
    """
    Print a skyplot with angle starting from North and going clockwise.
    
    Parameters
    ----------
    df: DataFrame
        The dataframe from the satellites
    Returns
    -------
    Skyplot 
    """
    # Only GLONASS
    if GLONASS==True and GPS==False:
        df = df.loc[df['system'] == 'GLONASS']
       
        az = df['azimuth'].to_list()
        el = df['elevation'].to_list()
    
    # Only GPS
    elif GPS==True and GLONASS==False:
        df = df.loc[df['system'] == 'GPS']
        
        az = df['azimuth'].to_list()
        el = df['elevation'].to_list()
    
    # Both
    elif GPS==True and GLONASS==True:
            
        az = df['azimuth'].to_list()
        el = df['elevation'].to_list()
        
    elif GPS==False and GLONASS==False:
        raise Exception("No constellation to plot")
        
    azrad=[np.pi/180*z for z in az]
    ax = plt.subplot(111, projection='polar')
    ax.set_ylim(bottom=90, top=0)
    ax.set_theta_zero_location("N")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise
    ax.set_title("Skyplot", va='bottom')
    
    ax.scatter(azrad,el)


    plt.show()

