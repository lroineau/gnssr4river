# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:55:23 2022

@author: souna
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle

def plotSnr(df, azrange=None, xaxis=0, show=False, save=True):
    """
    Subplot of SNR ratio to time. Possible to select only from satellites of 
    certain azimuth. If show=True, the subplot will be shown. Else it is saved
    in results directory. Note that due to large number of plot, it not adviced
    too set value to True.

    Parameters
    ----------
    df: DataFrame
        DataFrame from nmea.
    azrange: List
        window of azimuth
    xaxis: int
        1 for time and 0 for elevation

    Returns
    -------
    Plot

    """
    
    df = df.loc[((df.azsmth > azrange[0]) & (df.azsmth < azrange[1]))]

    row = int(np.ceil(len(df.index.get_level_values(2).unique())/3))

    fig, ax = plt.subplots(row, 3, figsize=(20, 10), dpi=80)

    l_PRN = df.index.get_level_values(2)

    # loop through the length of tickers and keep track of index
    n=0

    for val in enumerate(l_PRN.unique()):
        n+=1
        prn_val = val[1]
        
        # add a new subplot iteratively
        ax = plt.subplot(row, 3, n)

        if xaxis==1:
            df.xs(prn_val,level='PRN').reset_index().plot(x='time',y=['snr'],ax=ax)
        elif xaxis==0:
            df.xs(prn_val,level='PRN').reset_index().plot(x='elevsmth',y=['snr'],ax=ax)

        # just for legend, could be done directly with system...
        if prn_val<50:
            ax.legend([f"PRN {prn_val} GPS"])
        else:
            ax.legend([f"PRN {prn_val} Glonass"])
        ax.set_ylabel('SNR (dB/Hz)')
            
        if xaxis==1:
            fig.suptitle('SNR ratio to time for each satellite in file')
        elif xaxis==0:
            fig.suptitle('SNR ratio to elevation for each satellite in file')
            
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    if show==False and save==True:
        plt.close(fig)
        if xaxis==1:
            plt.savefig("SNR_elev.png")  
        elif xaxis==0:
            plt.savefig("SNR_time.png")  
        
    elif show==True and save==False:
        plt.show()
    
    elif show==True and save==True:
        plt.show()    
        if xaxis==1:
            plt.savefig("SNR_elev.png")  
        elif xaxis==0:
            plt.savefig("SNR_time.png")  
    
    
def detrende_signal(df,PRN,order=4):
    """
    Detrende SNR data using low-order polynomial.
    
    Parameters
    ----------
    df: DataFrame
        DataFrame from nmea.
    PRN: int
        PRN number of satellite to analyze
    order: int
        Order of polynomial trend to remove from the direct signal

    Returns
    -------
    dfp: DataFrame
        Rearanged DataFrame
    """

    # PRN number
    dfp = df.xs(PRN,level='PRN').reset_index()
    
    # Convert SNR to a linear scale (dB/Hz -> Volt/Volt) and calculate sin(elevation)
    dfp["sinelev"]=np.sin(dfp.elevsmth*np.pi/180)
    dfp['snrV'] = 10**(dfp.snr/20)

    # Detrend the signal
    p = np.poly1d( np.polyfit(dfp['sinelev'],dfp['snrV'],order))
    t = np.linspace(0, 1, len(dfp['snrV']))
    dfp['snrV'] = dfp['snrV']-p(t)

    # Remove data above 30°
    dfp = dfp.loc[dfp.sinelev < 0.4]
    
    return dfp
    
def plot_detrende(dfp,PRN):
    """
    plot detrende signal.
    
    Parameters
    ----------
    df: DataFrame
        DataFrame from nmea.
    PRN: int
        PRN number of satellite to analyze.
    Returns
    -------
    Plot
    """
    ax = plt.subplot()
    leg=[]
    ax.plot(dfp.sinelev,dfp.snrV)


    ax.set_ylabel('SNR [Volt/Volt]')
    ax.set_xlabel('sin(Elevation angle)')
    ax.set_title(f'SNR data with direct signal removed, PRN{PRN}') 

    return

def height_LSP(dfp,minH,maxH,PRN):
    """
    Use LSP to calculate the height.
    
    Parameters
    ----------
    dfp: DataFrame
        Rearanged DataFrame.
    minH, maxH: int
        Window of expected value for the height.
    PRN: int
        PRN number of satellite to analyze.
    Returns
    -------
    frequency,height: list
        list of value of LSP frequencies and their equivalent in height.
    maxF,maxAmp: float
        max value for frequency and amplitude of the LSP (i.e. the height).
    """
    # LSP
    frequency, power = LombScargle(dfp['sinelev'],dfp['snrV']).autopower()

    # Change x-axis to height instead of frequency
    cf = 299792458/1575.42e6
    height = frequency*cf/2 

    ij = np.argmax(power)
    maxF = height[ij]
    maxAmp = np.max(power)

    plt.plot(height, power) 
    plt.axvline(x=maxF, color='r', linestyle='-')
    plt.xlim(minH,maxH)
    plt.xlabel("Reflector height (meter)")
    plt.ylabel("Amplitude")
    plt.title("LSP for PRN{PRN}")
    
    return frequency,height,maxF,maxAmp
