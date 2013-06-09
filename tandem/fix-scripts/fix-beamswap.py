#!/usr/bin/env python

"""
beamswap.py
----------------

Swaps beam 11 and 13 around.

"""

import pyfits, numpy, pylab, sys, os, re
import pylab as plt

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
    print "Swapping beams 11 and 13..."
    sd_len = fits[1].data.shape[0]
    for row in range(sd_len):
        if fits[1].data["BEAM"][row] == 11:
            print ".",
	    fits[1].data["BEAM"][row] = 13
        elif fits[1].data["BEAM"][row] == 13:
            print "+",
            fits[1].data["BEAM"][row] = 11   
    fits.flush()
    fits.close()

