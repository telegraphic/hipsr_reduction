#!/usr/bin/env python

print "HIPSR PLOT GIBBS"
print "----------------"

import pyfits as pf, sys, os, re
import pylab as plt
import numpy as np

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


mb = sdFits('/home/dprice/data/P669/mbcorr/sdf/2012-10-29_1113-P669_g02_1285_p669_dec02a.sdfits')
hi = sdFits('/home/dprice/data/P669/hipsr/raw/hipsr_2012-10-29_1113-P669_g02_1285_p669_dec02a.sdfits')
snum = 32 # Sample numbers 

print "HIPSR:  %s"%hi.filename
print "MBCORR: %s"%mb.filename

hi_data_y    = hi.data[0,0,0,0,2900:3000]
mb_data_y    = mb.data[0,0,0,0,:]


plt.subplot(211)
plt.plot(mb.freqs, mb_data_y, color='#333333', label='MBCORR')
plt.xlim(1267.5,1269.0)
plt.ylim(33,60)
plt.ylabel("Uncalibrated Power [-]")
plt.legend(frameon=False)
plt.xticks([i for i in np.arange(1267.5,1269.0+1, 0.5)],[str(i) for i in np.arange(1267,1270+1, 0.5)])

plt.subplot(212)
plt.plot(hi.freqs[2900:3000]-2*0.0488, hi_data_y*15, color='#333333', label='HIPSR')
plt.xlim(1267.5,1269.0)
plt.ylim(28,60)
plt.ylabel("Uncalibrated Power [-]")
plt.xlabel("Frequency [MHz]")
plt.legend(frameon=False)
plt.xticks([i for i in np.arange(1267.5,1269.0+1, 0.5)],[str(i) for i in np.arange(1267,1270+1, 0.5)])

plt.show()


