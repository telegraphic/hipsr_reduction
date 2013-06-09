#!/usr/bin/env python

"""
beamswap.py
----------------

Swaps beam 11 and 13 around.

"""

import pyfits, numpy, pylab, sys, os, re
import pylab as plt
from optparse import OptionParser

# Basic option parsing 
p = OptionParser()
p.set_usage('fix-P669-all.py [filepath]')
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
        
        # Fix beam mapping (11 <->13)
        if fits[1].data["BEAM"][row] == 11:
            print ".",
	    fits[1].data["BEAM"][row] = 13
        elif fits[1].data["BEAM"][row] == 13:
            print "+",
            fits[1].data["BEAM"][row] = 11 

        # Reverse data if required
        if int(fits[1].data['CRVAL1'][0]) == int(1375e6):
            fits[1].data["DATA"][row,0,0,0] = fits[1].data["DATA"][row,0,0,0,::-1]
            fits[1].data["DATA"][row,0,0,1] = fits[1].data["DATA"][row,0,0,1,::-1]
        
        # Basic RFI flagging    
        flagged = fits[1].data["FLAGGED"][row,0,0,0,:] 
        flagged[0:1250]    = 1 # Out of band
        
        if int(fits[1].data['CRVAL1'][0]) == int(1325e6):
            #print "1325MHz"
            flagged[7400:]     = 1 # Strong RFI source and Thuraya to end of band
            flagged[6650:6690] = 1 # 1450MHz thing
            
            
        if int(fits[1].data['CRVAL1'][0]) == int(1375e6):
            #print "1375MHz"
            flagged[7000:]     = 1 # Strong RFI source and Thuraya to end of band
            flagged[5620:5660] = 1 # 1450MHz thing
        
        fits[1].data["FLAGGED"][row,0,0,0,:] = flagged
        fits[1].data["FLAGGED"][row,0,0,1,:] = flagged
    
    print "OK.\n"                       
    fits.flush()
    fits.close()

