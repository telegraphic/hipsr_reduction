#!/usr/bin/env python

"""
fix-P669-rfi-simple.py
----------------

Flag ubiquitous RFI

"""

import pyfits, numpy, pylab, sys, os, re, time
import pylab as plt
import numpy as np
path = os.getcwd()

# Regular expression to match SDFITS extension
regex = '([0-9A-Za-z-_]+).sdfits'
filelist = []

t = raw_input("Warning: this file overwrites data. Be careful. Press Enter.")

for filename in os.listdir(path):
    match = re.search(regex, filename)
    if match:
        filelist.append(filename)
    
for iname in filelist:
    name = os.path.join(path, iname)
    fits     = pyfits.open(name,mode='update')
    print "File: %s, c. freq: %s"%(iname, fits[1].data['CRVAL1'][0])
    print "Flagging known strong RFI sources...", 
    sd_len = fits[1].data.shape[0]
    freq = fits[1].data['CRVAL1'][0]
    if int(freq) == int(1325e6): print "1325MHz"
    if int(freq) == int(1375e6): print "1375MHz"
    for row in range(sd_len):
        
        flagged = fits[1].data["FLAGGED"][row,0,0,0,:] 
        flagged[0:1250]    = 1 # Out of band
        if int(freq) == int(1325e6):
            #print "1325MHz"
            flagged[7400:]     = 1 # Strong RFI source and Thuraya to end of band
            flagged[6650:6690] = 1 # 1450MHz thing
            
            
        if int(freq) == int(1375e6):
            #print "1375MHz"
            flagged[7000:]     = 1 # Strong RFI source and Thuraya to end of band
            flagged[5620:5660] = 1 # 1450MHz thing
        
        fits[1].data["FLAGGED"][row,0,0,0,:] = flagged
        fits[1].data["FLAGGED"][row,0,0,1,:] = flagged
    
    fits.flush()
    fits.close()

