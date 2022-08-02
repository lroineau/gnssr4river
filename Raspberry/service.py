#!/usr/bin/env python
"""
This file is to use in case you want to write data from the antenna during a 
survey without access to ethernet. It simply open a new file every time the 
raspberry boots. See how to built service file for more information:
link_to_git

Author: Lubin Roineau
"""

import serial
import os
import gzip
import yaml

# Read content of config file 
with open("config.yml", "r") as yamlfile:
    data = yaml.load(yamlfile, Loader=yaml.FullLoader)

port = data[2]
path_to_file = data[3]
name = data[4]

# Serial port
ser = serial.Serial('/dev/ttyAMA0', baudrate=9600)    
    
# Setting up increment to check file name
c = 0
    
# Check if there is already a file
while os.path.exists(path_to_file+ '/{}_{}.gz'.format(name,c)) is True:
    f = gzip.open(path_to_file+ '/{}_{}.gz'.format(name,c),'r')
    f.close()
    c+=1

# Create file
filenamegz = path_to_file+ '/{}_{}.gz'.format(name,c)
f = gzip.open(filenamegz, 'wt')

# While raspberry on
while True:
    
    try:
        line = ser.readline()
        f.write(line.strip() + "\n")
        
    except KeyboardInterrupt:
        f.close()      
        break