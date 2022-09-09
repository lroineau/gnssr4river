#!/usr/bin/env python3
"""
This file is to use in case you want to write data from the antenna during a 
survey using a service file. For more information:
link_to_git

Author: Lubin Roineau
Modified: Roelof Rietbroek
"""

import serial
import os
import gzip
import yaml
from datetime import datetime
import re

def readconfig():
    """Read the content of configuration file"""
    configfile=os.path.join(os.path.expanduser('~'),"nmeaconfig.yml")
    with open(configfile, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)
    return cfg


def getlogstream(filebase):

    # Setting up increment to check file name
    c = 0
    
    # Check if there is already a file
    filenamegz=f"{filebase}_{c:02d}.gz"
    while os.path.exists(filenamegz) is True and c < 99:
        c+=1
        
        # Create file name
        filenamegz=f"{filebase}_{c:02d}.gz"
        
        
    if c <=99:
        print("opening",filenamegz) 
        return gzip.open(filenamegz, 'wb')
    else:
        #limit reached
        return None

def getdate(line):
    return datetime.strptime(line.split(b",")[9].decode('utf-8'),"%d%m%y")

def main():
    conf=readconfig()
    fid=None
    prevdate=datetime.min.date()
    #regular expression matching an RMC message
    rmcregex=re.compile(b'^\$G[NPL]RMC')
    
    #open serial port
    ser = serial.Serial(conf['device'],baudrate=conf['baudrate'])
    while True:
        try:
            line=ser.readline()

            if rmcregex.match(line):
                currentdate=getdate(line)
                #possibly open a new logger stream on date changes
                if prevdate < currentdate.date():
                    if fid:
                        fid.close()
                    filenamebase=os.path.join(conf['data_dir'],f"{conf['file_base']}_{currentdate.date().isoformat()}")
                    fid=getlogstream(filenamebase)
                    prevdate=currentdate.date()

            if fid:
                fid.write(line)
            else:
                #continue until a proper date has been extracted from the RMC message
                continue

        except KeyboardInterrupt:
            if fid:
                fid.close()


if __name__ == "__main__":
    main()
