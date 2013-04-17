#!/usr/bin/env python

print "HIPSR BANDPASS FITTER"
print "---------------------"

import pyfits as pf, sys, os, re
import pylab as plt
import numpy as np
import pywt
from scipy import interpolate
mode = pywt.MODES.sp1

from numpy import log10, abs, arange, array


def waveletExcise(xydata, threshold):
    """Excise RFI from a data using a discrete wavelet transform (DWT).
    A DWT is applied to the ydata, resulting in an approximated
    spectrum, cA, of len(xydata)/2, and detail coefficients, cD,
    of the same length. This is essentially lowpassed and highpassed
    data. High values in the detail coefficient vector cD correspond to
    large deviations in the spectrum - which we are assuming is unwanted.
    Any value in cD above the threshold will then be flagged as RFI, and
    removed from the dataset. 
    """
    
    xdata = xydata[:,0]
    ydata = xydata[:,1]
    
    # Do a discrete wavelet transformation on our data
    (cA, cD) = pywt.dwt(ydata, 'haar', mode='cpd')

    # Create flags for details coefficients over threshold value
    flags = []

    for i in range(0,len(cD)):
        if abs(cD[i]) > threshold:
            flags.append(i)
    flags.reverse()
    
    # Unfortunately, we can't delete numpy array elements by index
    # So we have to convert it to a list and then back again
    
    newxydata = []
    for i in range(0,len(xdata)):
        if(i % 2):
            row = [xdata[i],cA[(i-1)/2]]
            newxydata.append(row)
    
    for flag in flags: 
        del(newxydata[flag])
    
    # Now turn it back into an array and return
    newxydata = array(newxydata)
    
    numflagged = len(flags)
    numdatapoints = len(newxydata)
    print "Number of flagged points: %i"%numflagged
    print "Number of remaining datapoints: %i"%numdatapoints
    
    return newxydata, flags


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


hi = sdFits('/home/dprice/data/P669/hipsr/raw/hipsr_2012-10-29_1113-P669_g02_1285_p669_dec02a.sdfits')

snum = 32 # Sample numbers 

print "HIPSR:  %s"%hi.filename

data_x    = hi.data[0,0,0,0,:]
freqs     = hi.freqs

xydata = np.column_stack((freqs, data_x))
xydata, flags = waveletExcise(xydata, 0.005)

tck = interpolate.splrep(xydata[:,0],xydata[:,1],s=0)
data_s = interpolate.splev(freqs,tck,der=0)

plt.subplot(311)
plt.plot(freqs, data_x, color='#333333', label='Original data')
plt.xlim(1150,1500)
plt.ylim(-0.09,6)
plt.legend()

plt.subplot(312)
plt.plot(xydata[:,0], xydata[:,1], color='#333333', lw=1.5, label='Data post DWT flagging')
plt.plot(freqs, data_s, color='#CC0000', label='Interpolated spline')
plt.xlim(1150,1500)
plt.ylim(-0.09,6)
plt.legend()

plt.subplot(313)
plt.plot(freqs, data_x/data_s, color='#333333', label='Final bandpass')
plt.xlim(1150,1500)
plt.ylim(0.645,1)
plt.xlabel("Frequency [MHz]")
plt.legend()

plt.show()

plt.plot(freqs, data_x/data_s, color='#333333', label='Final bandpass')
plt.xlim(1415,1425)
plt.xlabel("Frequency [MHz]")
plt.ylim(0.7,0.85)
plt.legend()
plt.show()

