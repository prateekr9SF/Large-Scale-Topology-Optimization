#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 00:36:41 2025

@author: wz10
"""

import numpy as np
import sys
import os
import subprocess
def runcalTop(x,File, penalty, rmin):
    with open("density.dat", mode="w+") as file:
        for xi in x:
            file.write(str(xi)+"\n")
    File=File
    cmd = ["calTop.exe", File, "-p", str(penalty), "-r", str(rmin)]
    result = subprocess.run(cmd)
    
#Add stuff to read penal rmin and x 
config = sys.argv[1]
#config = 'OPTIM/DIRECT/config_tmpl.txt'
fid = open(config,"r")
lines = fid.readlines()
fid.close()

data = lines[0][0:-1]
x = np.fromstring(lines[1], sep=',')
fid = open(data,"r")
lines = fid.readlines()
fid.close()

nDV = int(lines[0])
volfrac = float(lines[1])
rmin = float(lines[2])
penalty = int(lines[3])
File = lines[4].rstrip("\n")
NCPU = int(lines[5])
#Write the density file
runcalTop(x, File, penalty, rmin)

#Run the actual solver

    
