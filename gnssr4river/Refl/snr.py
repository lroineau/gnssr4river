# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:55:23 2022

@author: souna
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def show(df):
    dflow=df.loc[(df.elevsmth < 40) & (df.elevsmth > 2) &((df.azsmth > 200) & (df.azsmth < 230))]
    # dflow=df.where((df.elevsmth < 25) & (df.elevsmth > 4) & (df.snr.notna()) )
    dflow["sineelev"]=np.sin(dflow.elevsmth*np.pi/180)
    i=0
    ax=plt.subplot()
    leg=[]
    for name,grp in dflow.groupby(["system","PRN","segment"]):
        ax.plot(grp.elevsmth,grp.snr)
        leg.append(f"{name[0]} PRN{name[1]}")
        
        ax.set_ylabel('Carrier to noise [dBhz]')
        ax.set_xlabel('Elevation angle [degrees]')
        ax.legend(leg)
        ax.set_title('Omleidingskanaal reflections from 200 < azimuth < 230') 
        

def CNR2SNR(df):
    """
    Takes the Carier to Noise Ratio from the nmea file and return the SNR.

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    dflow["sineelev"]=np.sin(dflow.elevsmth*np.pi/180)
    dflow.groupby(["system","PRN","segment"])
    grp.elevsmth,grp.snr
    
    SNR = 10**(CNR/20)
    ax=plt.subplot()
    
    for name,grp in dflow.groupby(["system","PRN","segment"]):
        ax.plot(grp.sineelev,grp.snr)        
        ax.set_ylabel('SNR [Volt/Volt]')
        ax.set_xlabel('sin(elevation) [degrees]')
        ax.set_title('Omleidingskanaal reflections from 200 < azimuth < 230') 
    