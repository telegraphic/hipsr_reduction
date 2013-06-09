#!/usr/bin/env python

print "HIPSR/MBCORR data checker"
print "-------------------------"

print """
This script computes & compares the RMS statistics of HIPSR and MBCORR data:

* Two matching SD-FITS files are opened (the HIPSR and MBCORR)
* The central 64 channels of each spectrum are extracted
* The RMS levels of these 64 channels are computed
* The average RMS level of all spectra in the file are computed (RMS AVG)
* The stdev is also computed (RMS STD)
* The ratio of HIPSR:MBCORR RMS is then shown for each channel (it should be ~1)

"""

import pyfits as pf, sys, os, re
import pylab as plt
import numpy as np


def avgStd(data):
    data = np.array(data)
    avg, std  = np.average(data), np.std(data)
    return avg, std

class sdFits(object):
    def __init__(self, filename):
        fits     = pf.open(filename)
       
        self.filename     = fits.filename()
        self.time    = fits[1].data['TIME']
        self.date    = fits[1].data['DATE-OBS']
        self.data    = fits[1].data['DATA']
        self.beam    = fits[1].data['BEAM']
        self.flagged = fits[1].data['FLAGGED']
        self.header  = fits[1].header
        self.shape   = np.shape(self.data)

        ref_pix   = fits[1].data['CRPIX1']  # Reference pixel
        ref_val   = fits[1].data['CRVAL1']  # Value at reference pixel (in Hz)
        ref_delt  = fits[1].data['CDELT1']  # Delta between pixels
        num_pix   = self.data.shape[-1]# Last axis of data array (multidimensional)
        
        # Might be better to look this up than assume it's col 25
        self.data_name = fits[1].header.get('TTYPE25')   
        self.data_unit = fits[1].header.get('TUNIT25')
        
        self.freqs = (np.arange(0,num_pix,1) * ref_delt[0] + ( ref_val[0] - ref_pix[0] * ref_delt[0] ) ) / 1e6 
        self.freq_type = fits[1].header.get('CTYPE1')
        self.freq_unit = 'MHz' # TODO: check this
        
        fits.close()

# Config
#mb = sdFits('/home/dprice/data/P641/SDFDATA_CAL/20121025/2012-10-25_1911-P641_west2_1315_P641.sdfits')
#hi = sdFits('/home/dprice/data/P641/hipsr_cal_sdfits/west/hipsr_2012-10-25_1911-P641_west2_1315_P641.sdfits')

mb = sdFits('/home/dprice/data/P641/SDFDATA_CAL/20121025/2012-10-25_1027-P641_east1_1315_P641.sdfits')
hi = sdFits('/home/dprice/data/P641/hipsr_2012-10-25_1027-P641_east1_1315_P641.sdfits')
snum = 32 # Sample numbers 

print "HIPSR:  %s"%hi.filename
print "MBCORR: %s"%mb.filename

print "\n"    
print "      |       HIPSR       |       MBCORR      |         "
print "POL A | RMS AVG   RMS STD | RMS AVG   RMS STD | RATIO   "
print "--------------------------------------------------------"

for beamcnt in range(1,13+1):
    ratios, hi_rms, mb_rms = [], [], []
    for row in range(hi.shape[0]):   
        if hi.beam[row] == beamcnt:
            hi_data_x    = hi.data[row,0,0,0,:]
            hi_flagged_x = hi.flagged[row,0,0,0,:]

            mb_data_x    = mb.data[row,0,0,0,:]
            mb_flagged_x = mb.data[row,0,0,0,:]
            
            hi_std = np.std(hi_data_x[len(hi_data_x)/2-snum:len(hi_data_x)/2+snum])        
            mb_std = np.std(mb_data_x[len(mb_data_x)/2-snum:len(mb_data_x)/2+snum])
            
            hi_rms.append(hi_std)
            mb_rms.append(mb_std)
  
    hi_avg, hi_std = avgStd(hi_rms)
    mb_avg, mb_std = avgStd(mb_rms)
    
    if all((mb_avg > 1e-9, hi_avg > 1e-9)):
        ratio = hi_avg / mb_avg
    else:
        ratio = 0
    
    print "%5i | %2.4f +/- %2.4f | %2.4f +/- %2.4f | %2.4f"%(beamcnt, hi_avg, hi_std, mb_avg, mb_std, ratio)

print "\n"    
print "      |       HIPSR       |       MBCORR      |         "
print "POL B | RMS AVG   RMS STD | RMS AVG   RMS STD | RATIO   "
print "--------------------------------------------------------"

for beamcnt in range(1,13+1):
    ratios, hi_rms, mb_rms = [], [], []
    for row in range(hi.shape[0]):   
        if hi.beam[row] == beamcnt:
            hi_data_y    = hi.data[row,0,0,1,:]
            hi_flagged_y = hi.flagged[row,0,0,1,:]

            mb_data_y    = mb.data[row,0,0,1,:]
            mb_flagged_y = mb.data[row,0,0,1,:]
 
            hi_std = np.std(hi_data_y[len(hi_data_y)/2-snum:len(hi_data_y)/2+snum])        
            mb_std = np.std(mb_data_y[len(mb_data_y)/2-snum:len(mb_data_y)/2+snum])
            
            hi_rms.append(hi_std)
            mb_rms.append(mb_std)

    
    hi_avg, hi_std = avgStd(hi_rms)
    mb_avg, mb_std = avgStd(mb_rms)
    
    if all((mb_avg > 1e-9, hi_avg > 1e-9)):
        ratio = hi_avg / mb_avg
    else:
        ratio = 0
    
    print "%5i | %2.4f +/- %2.4f | %2.4f +/- %2.4f | %2.4f"%(beamcnt, hi_avg, hi_std, mb_avg, mb_std, ratio)   
    
