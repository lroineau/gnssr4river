# -*- coding: utf-8 -*-
"""
Read nmea file. This code is part of frommle2 made by Roelof Rietbrock. See more here:
https://github.com/strawpants/frommle2

@author: Lubin Roineau, ENSG-Geomatics (internship at UT-ITC Enschede), Aug 26, 2022
"""

# Imports
import gzip as gz
from datetime import datetime,timedelta
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

###############################################################################

cconv=(1-100/60)
swopfac={"E":+1,"W":-1,"N":1,"S":-1}
def parseDeg(deg,EWNS):
    decdeg=deg/60+cconv*int(deg/100)
    return swopfac[EWNS]*decdeg

def smoothDegrees(degarray,timev,irec=0):
        """ Smooths degree array which only have degree resolution to a version which varies more smoothly (i.e. no jumps)"""
        
        if len(degarray) == 1:
            #Corner case just return original value
            return degarray


        ddif=np.diff(degarray)

        if np.count_nonzero(ddif) == 0:
            #no need to interpolate in this case as all values are the same
            return degarray

        
        if irec > 1:
            assert("This function is not supposed to be called with a recursion depth more than 1)")
            #Possibly split into contigous sections (i.e. crossing  0-360 border or jumping in time)
        dsections=np.argwhere(np.logical_or((np.abs(ddif) > 180.0), (np.diff(timev) > np.timedelta64(30,'s'))))  
        if dsections.size != 0:
            degsmth=np.full([len(degarray)],np.nan)
            for ist,iend in zip(np.insert(dsections+1,0,0),np.append(dsections+1,len(degarray))):
                #call this function for the sections (it shouldn;t execute this if part)
                degsmth[ist:iend]=smoothDegrees(degarray[ist:iend],timev[ist:iend],irec=irec+1)
        else:
            
            #compute for the contiguous section
            #detect the jumps (and add the first and final index) value
            ijumps=np.insert(np.append(np.argwhere(ddif != 0).squeeze()+1,len(degarray)),0,0)

            #what are half the interval lengths between the jump locations (we need to shift the interpolant input by this)
            iintHalfDistance=np.diff(ijumps)/2
            #create the support points for the interpolant
            isup=ijumps[0:-1]+iintHalfDistance
            degsup=degarray[ijumps[0:-1]]
            inonnan=~np.isnan(degarray[ijumps[0:-1]])
            #create an interpolant based on the values where a jump was detected
            iinterp=interp1d(isup[inonnan],degsup[inonnan],kind='linear',fill_value="extrapolate")
            degsmth=iinterp([float(i) for i in range(len(degarray))])
            

            degsmth[np.isnan(degsmth)]=np.nan #refill original nan values with nans again
        
        return degsmth

def parseGNRMC(ln):
    
    spl=ln.split(",")
    if spl[2] == "V":
        return {}
    hr=int(spl[1][0:2])
    mn=int(spl[1][2:4])
    sec=float(spl[1][4:])
    date=datetime.strptime(spl[9],"%d%m%y")
    dt={"time":date+timedelta(hours=hr,minutes=mn,seconds=sec)}
    #parse weird DDDMM.MMMMM format
    dt["lat"]=parseDeg(float(spl[3]),spl[4])
    dt["lon"]=parseDeg(float(spl[5]),spl[6])
    return dt

def parseGNGSV(ln):
    #split line without considering the checksum at the end
    spl=ln[0:-4].split(",")
    dt={}
    system=spl[0][1:3].replace('GL','GLONASS').replace('GP','GPS')
    #loop over available satellite prn's
    for i in range(4,len(spl),4):
        try:
            prn=f"PRN{spl[i]}"
            elev=float(spl[i+1])
            az=float(spl[i+2])
            snr=float(spl[i+3])
        except ValueError:
            #It may be possible that ,,, entries occur, so we'll just ignore those
            continue

        dt[prn]={"system":system,"elev":elev,"az":az,"snr":snr}
    
    return dt   


dispatchParse={"$GPRMC":parseGNRMC,"$GPGSV":parseGNGSV,"$GNRMC":parseGNRMC,"$GLGSV":parseGNGSV}

###############################################################################

def readnmea(fidorfile):
    """Parses a nmea file/stream and puts the output in a pandas dataframe"""
    if type(fidorfile) == str:
        if fidorfile.endswith('.gz'):
            fid=gz.open(fidorfile,'rt')
        else:
            fid=open(fidorfile,'rt')
    else:
        fid=fidorfile



    #loop over the buffer and parse messages as we go
    nmeacycle={}
    nmeadata=[]
    for ln in fid:
        if ln.startswith("$"):
            try:
                nmeacycle.update(dispatchParse[ln[0:6]](ln))

                #possibly append this cycle data to nmeadata when a time tag is present
                if "time" in nmeacycle and (sum(k.startswith("PRN") for k in nmeacycle.keys()) > 0):

                    basedict={k:v for k,v in nmeacycle.items() if not k.startswith("PRN")}

                    #unwrap the different PRN's into separate rows
                    for ky,val in nmeacycle.items():
                        if ky.startswith("PRN"):
                            nmeadata.append({**basedict,"PRN":int(ky[3:]),**val})
                    #reset nmeacycle dict
                    nmeacycle={}
            except KeyError:
                continue

    #create a dataframe and set multiindex
    df=pd.DataFrame(nmeadata)
    #We wan to get rid of rows which don;t have a elevation of azimuth in them
    # df.dropna(subset=["elev","az"],inplace=True)
    

    #since the resolution of the elevatio, and azimuthn is only on discrete degrees, let's create smoothed versions (better for plotting etc)
    # Furthermore, we also want to identify the different ascending/descending segment, even/oddly numbered repspectively so we can select of those too
    df["elevsmth"]=np.nan 
    df["azsmth"]=np.nan 
    df["segment"]=-1
    for name,grp in df.groupby("PRN"):
        time=grp.time.to_numpy()
        for dvi,dvo in [("elev","elevsmth"),("az","azsmth")]:
            dsmooth=smoothDegrees(grp[dvi].to_numpy(),time)
            #put the stuff back in the dataframe
            df.loc[grp.index,dvo]=dsmooth


            if dvi == "elev" and len(dsmooth) > 2:
                #also figure out the valid iascending (will be even numbered, descending will be oddly numbered) segments per PRN
                pdif=(np.diff(dsmooth)>0).astype(float)
                segment=np.insert(np.cumsum(np.abs(np.diff(pdif))),0,[0,0]).astype(int)
                if pdif[0] == 0.0:
                    #starts with a descending node
                    #make sure it's oddly numbered
                    segment=segment+1

                df.loc[grp.index,"segment"]=segment

                df=df.loc[df.snr.notna()]
                
    return df.set_index(["time","system","PRN","segment"])