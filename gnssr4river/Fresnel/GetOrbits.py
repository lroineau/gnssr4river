# -*- coding: utf-8 -*-
"""
This code is to retrieve the orbits of GPS and GLONASS satellites

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import wget
import os
import numpy as np
import unlzw3
from pathlib import Path
import subprocess
import pandas as pd
import datetime as d

###############################################################################

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

###############################################################################

def retrieve_orbits(year=None, month=None, day=None, hour=None, minute=None, second=None):
    """
    Retrieve the orbits of GPS and Glonass constellations from ESA site. If date is not specified, takes current date.

    Parameters
    ----------
    year, month, day, hour, minute, second : int
        Year (xxx), month, day, hour, minute, second (default is current daytime)
    
    Return
    ------
    
    """
    # Calcul the current date
    cur_date = d.datetime.now()
    cur_date = cur_date.strftime("%Y,%m,%d,%H,%M,%S")
    cur_date = cur_date.split(",")
    cur_year, cur_month, cur_day, cur_hour, cur_minute, cur_second = int(cur_date[0]), int(cur_date[1]), int(cur_date[2]), int(cur_date[3]), int(cur_date[4]), int(cur_date[5])

    # if the date is not given by the user, take current date
    if year==None:
        year, month, day, hour, minute, second = cur_year, cur_month, cur_day, cur_hour, cur_minute, cur_second
        
    # if user gives an incorrect date
    giv_date = str(year + ',' + month + ',' + day + ',' + hour + ',' + minute + ',' + second)
    if giv_date > cur_date:
        raise Exception("Cannot give a date that is yet to come !")  

    # Retrieve first the gps week
    GPS_wk, GPS_sec_wk = gpsweek(year, month, day, hour, minute, second)
    
    # Create directory
    dirName = 'Orbits'
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory {} created".format(dirName))
    else:    
        print("Directory {} already exists".format(dirName))
        
    # Name of file to retrieve
    # Zipped name
    filenameZ = 'esu{}'.format(GPS_wk) + str(int(GPS_sec_wk/86400)) + '_00.sp3.Z' 
    # Unzipped name    
    filename = 'esu{}'.format(GPS_wk) + str(int(GPS_sec_wk/86400)) + '_00.sp3'         
    # data link   
    url = 'http://navigation-office.esa.int/products/gnss-products/{}/{}'.format(GPS_wk,filenameZ)
    # Check if file already exists 
    if os.path.exists('Orbits/{}'.format(filename)) is True:
        print("File already exists")
    # Else it is downloaded
    else:
        wget.download(url, out=dirName)
        subprocess.call(['uncompress',filename])
        subprocess.call(['mv',filename, dirName])
        if os.path.exists('Orbits/{}'.format(filename)) is True:
            print("Data download success")
        else:
            print("Fail to retrieve data")

    # uncompress file
    # subprocess.call(['uncompress', 'Orbits/esu{}0_00.sp3.Z'.format(GPS_wk)])
    # orb = unlzw3.unlzw(Path('Orbits/esu{}0_00.sp3.Z'.format(GPS_wk)))

    return 

###############################################################################  

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