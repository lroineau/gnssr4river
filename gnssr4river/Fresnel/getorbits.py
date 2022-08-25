# -*- coding: utf-8 -*-
"""
This code is to retrieve the orbits of GPS and GLONASS satellites
@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

### Imports ###
import requests
import os
import numpy as np
import unlzw3
from pathlib import Path
import pandas as pd
from datetime import datetime 
from io import BytesIO
import geod
import matplotlib.pyplot as plt

###############################################################################

def gpsweek(date=datetime.now()):
    """
    Return the gps week and secondes for a given date. 
    
    Parameters
    ----------
    date: datetime
        Date used for determining the gps week.
    
    Return
    ------
    GPS_wk, GPS_sec_wk : int       
        The gps week and second of the week.
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
    orb: bytes
        sp3 file.
    """
    # if user gives an incorrect date
    if date > datetime.now():
        raise RuntimeError("Cannot give a date that is yet to come!")  

    # Retrieve first the gps week
    GPS_wk, GPS_sec_wk = gpsweek(date)    
    t = str(int(GPS_sec_wk/86400))

    # Create directory
    dirName = 'Orbits'
    # Create target Directory if don't exist
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print(f"Files will be stored in {dirName}")
        
    # List of possible file updated
    update=['18','12','06','00']
    a=0
    for i in range(len(update)):
        # Try to see which one is the latest file
        xx=update[a]
        filenameZ = f'esu{GPS_wk}' + t + '_'+xx+'.sp3.Z' 
        # If it exists localy, no need to download it
        if os.path.exists('Orbits/{}'.format(filenameZ)) is True:
            print('File already exists')
            orb = unlzw3.unlzw(Path(f'Orbits/{filenameZ}'))
            return orb
        # If not, download it
        else:
            # data link   
            url = f'http://navigation-office.esa.int/products/gnss-products/{GPS_wk}/{filenameZ}'
            r = requests.get(url)
            # Check the latest file updated for the given date
            if r.status_code == 200:
                open(dirName +'/' +filenameZ, 'wb').write(r.content)
                # Data download success
                if os.path.exists(f'Orbits/{filenameZ}') is True:
                    print("Data download success")
                    orb = unlzw3.unlzw(Path(f'Orbits/{filenameZ}'))
                    return orb
                # Data download failled
                if os.path.exists(f'Orbits/{filenameZ}') is False:
                    print("Fail to retrieve data") 
            else:
                a+=1

###############################################################################  

def read_sp3(file,predict_only=False):
    """
    Read the sp3 file and turn it to a DataFrame. Can also be a Bytes file if unzipped directly with unlzw3.
    
    Parameters
    ----------
    file: String or Bytes
        The sp3 file or the local unzipped file in Python.
    predict_only: Boolean
        Set to True if only predicted orbits are wanted.
        
    Return
    ------
    df_sp3: DataFrame
        The content of the sp3 on a DataFrame.
    """
    # See if the file is a sp3 on hard drive or directly a Bytes in Python
    if type(file)==bytes:        
        f = BytesIO(file) 
    else:
        f = open(f'Orbits/{file}')
    
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
        date = []
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
                date.append(dtepoch)
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
    df_sp3 = pd.DataFrame({"date":date,
                           "system":sys,
                           "week":week.astype(int),
                           "tow":tow, 
                           "x":x,
                           "y":y,
                           "z":z,
                           "prn":prn.astype(int),
                           "clock":clock})
        
    # If only predicted orbits are kept
    if predict_only is True:
        df_sp3 = df_sp3.loc[(df_sp3['date'] >= datetime.now())]
    
    return df_sp3

###############################################################################

def clean_dir(dir):
    """
    Simple function to clean a directory.
    """
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))
    return

###############################################################################

def skyplot(df,GPS=True,GLONASS=True):
    """
    Print a skyplot with angle starting from North and going clockwise. 
    By default plot GLONASS and GPS but can be filtered by setting to False.
    
    Parameters
    ----------
    df: DataFrame
        The dataframe from the satellites.
    
    Returns
    -------
    Skyplot 
    """
    ax = plt.subplot(111, projection='polar')
    ax.set_ylim(bottom=90, top=0)
    ax.set_theta_zero_location("N")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise
    
    # Only GLONASS
    if GLONASS==True and GPS==False:
        df = df.loc[df['system'] == 'GLONASS']
       
        azgl = df['azimuth'].to_list()
        elgl = df['elevation'].to_list()
        
        azradgl=[np.pi/180*z for z in azgl]

        ax.scatter(azradgl,elgl,s=2,color='green')
        ax.set_title("Skyplot of GLONASS satellites (2 days)", va='bottom')

    
    # Only GPS
    elif GPS==True and GLONASS==False:
        df = df.loc[df['system'] == 'GPS']
        
        azgp = df['azimuth'].to_list()
        elgp = df['elevation'].to_list()
        
        azradgp=[np.pi/180*z for z in azgp]

        ax.scatter(azradgp,elgp,s=2,color='red')
        ax.set_title("Skyplot of GPS satellites (2 days)", va='bottom')

    # Both
    elif GPS==True and GLONASS==True:
        dfgp = df.loc[df['system'] == 'GPS']
        dfgl = df.loc[df['system'] == 'GLONASS']
        
        azgp = dfgp['azimuth'].to_list()
        elgp = dfgp['elevation'].to_list()
        azgl = dfgl['azimuth'].to_list()
        elgl = dfgl['elevation'].to_list()
        
        azradgl=[np.pi/180*z for z in azgl]
        azradgp=[np.pi/180*z for z in azgp]
        
        ax.scatter(azradgp,elgp,s=2,color='red',label='GPS')
        ax.scatter(azradgl,elgl,s=2,color='green',label='GLONASS')
        ax.legend(loc='lower left', bbox_to_anchor=(-0.2, -0.2),fancybox=True, shadow=True)
        ax.set_title("Skyplot of GLONASS and GPS satellites (2 days)", va='bottom')


    elif GPS==False and GLONASS==False:
        raise Exception("No constellation to plot")

    plt.show()
             