#!/usr/bin/env python

"""
beamswap.py
----------------

Swaps beam 11 and 13 around.

"""

import pyfits, sys, os, re
import pylab as plt
import numpy as np
from optparse import OptionParser

# Basic option parsing 
p = OptionParser()
p.set_usage('rfi-flagger.py [filepath]')
p.set_description(__doc__)
(options, args) = p.parse_args()

path = args[0]

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
    sd_len = fits[1].data.shape[0]
    freq = fits[1].data['CRVAL1'][0]

    # Change central frequency
    if int(freq) == int(1315e6):
        print "updating central frequency..."
        fits[1].data['CRVAL1'][:] = 1325e6
    
    # Main loop
    for row in range(sd_len):
        
        data_x    = fits[1].data['DATA'][row,0,0,0,:]
        data_y    = fits[1].data['DATA'][row,0,0,1,:]
        flagged_x = fits[1].data['FLAGGED'][row,0,0,0,:]
        flagged_y = fits[1].data['FLAGGED'][row,0,0,1,:]        
        
        # RFI flagging
        dx, dy, dxy = np.abs(data_x), np.abs(data_y), np.abs(data_x - data_y)
        sigma = 3.5 # Flagging threshold value
        for d in (dx,dy,dxy):
            dm = d[len(d)/2-len(d)/8:len(d)/2+len(d)/8] #assuming less RFI @ band centre
            threshold = np.average(dm) + sigma * np.std(dm)
            flagged_x[d > threshold] = 1
            flagged_y[d > threshold] = 1
        
        fits[1].data["FLAGGED"][row,0,0,0,:] = flagged_x
        fits[1].data["FLAGGED"][row,0,0,1,:] = flagged_y
    
    print "OK.\n"                       
    fits.flush()
    fits.close()

