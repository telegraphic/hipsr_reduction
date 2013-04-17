#!/usr/bin/env python

"""
sdfits_viewer.py
----------------

Viewer for Parkes multibeam SD-FITS files.

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
    imname = 'tmp'
        
    
    header   = fits[1].header
    
    # Frequency-like axis values
    ref_pix   = fits[1].data['CRPIX1']  # Reference pixel
    ref_val   = fits[1].data['CRVAL1']  # Value at reference pixel (in Hz)
    ref_delt  = fits[1].data['CDELT1']  # Delta between pixels
    num_pix   = 8192       # Last axis of data array (multidimensional)
    
    data_name = fits[1].header.get('TTYPE25')   # Might be better to look this up than assume it's col 25
    data_unit = fits[1].header.get('TUNIT25')
    
    freq_type = fits[1].header.get('CTYPE1')
    freq_unit = 'MHz'                           # Might be better to check this
    
    freqs = (numpy.arange(0,num_pix,1) * ref_delt[0] + ( ref_val[0] - ref_pix[0] * ref_delt[0] ) ) / 1e6 
    
    
    if int(fits[1].data['CRVAL1'][0]) == int(1315e6):
        print "updating central frequency..."
        fits[1].data['CRVAL1'][:] = 1325e6
        fits.flush()

    if int(fits[1].data['CRVAL1'][0]) == int(1375e6):
        print "Reversing data order..."

        sd_len = fits[1].data.shape[0]
        for row in range(sd_len):
            fits[1].data["DATA"][row,0,0,0] = fits[1].data["DATA"][row,0,0,0,::-1]
            fits[1].data["DATA"][row,0,0,1] = fits[1].data["DATA"][row,0,0,1,::-1]
        fits.flush()
    fits.close()

