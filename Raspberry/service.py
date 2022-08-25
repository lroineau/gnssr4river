#!/usr/bin/env python
"""
This file is to use in case you want to write data from the antenna during a 
survey using a service file. For more information:
link_to_git

Author: Lubin Roineau
"""

import serial
import os
import gzip
import yaml
import datetime

# Read content of config file 
with open("config.yml", "r") as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
    
    hour= cfg['service']['daily']
    port = cfg['parameters']['port']
    path_to_file = cfg['parameters']['path_to_file']
    name = cfg['parameters']['file_name']


# Serial port
ser = serial.Serial('/dev/ttyAMA0', baudrate=port)    


###############################################################################
#                               For daily files:                              #
############################################################################### 
if hour is True:
    
    def checkfile():
        now = datetime.now()
        date = now.strftime('%Y-%m-%d')

        filenamegz = '{}_{}.gz'.format(name,date)
        
        # Check if the file already exists
        if os.path.exists(path_to_file+'/'+filenamegz) is True:
            # If yes, continue to write in it
            f = gzip.open(filenamegz, 'wt')
        # Else, open a new one
        else:
            filenamegz = path_to_file+'/'+filenamegz
            f = gzip.open(filenamegz, 'wt')
        return f
    
    while True:
    
        try:
            
            f = checkfile()
                
            while datetime.now().strftime('%H:%M') != '00:00':

                line = ser.readline()
                f.write(line.strip() + "\n")
        
            f.close()
            
        except KeyboardInterrupt:
            f.close()      
            break             
        
###############################################################################
#                      For simple acquisition with increment                  #
###############################################################################           
elif hour is False:
    
    # Setting up increment to check file name
    c = 0
    
    # Check if there is already a file
    while os.path.exists(path_to_file+ '/{}_{}.gz'.format(name,c)) is True:
        c+=1
        
        # Create file
        filenamegz = path_to_file+ '/{}_{}.gz'.format(name,c)
        f = gzip.open(filenamegz, 'wt')
        
    while True:
    
        try:
            line = ser.readline()
            f.write(line.strip() + "\n")
        
        except KeyboardInterrupt:
            f.close()      
            break