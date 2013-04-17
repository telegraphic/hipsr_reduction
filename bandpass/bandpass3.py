#!/usr/bin/env python

print "HIPSR BANDPASS FITTER"
print "---------------------"

import pyfits as pf, sys, os, re, time
import pylab as plt
import numpy as np
import pywt
from scipy import interpolate
from scipy.ndimage.filters import convolve1d

from numpy import log10, abs, arange, array

def timer(func):
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper

def waveletDenoise(data, threshold=1, wavelet='haar'):
    """ Wavelet denoise test """
    # Decompose
    wl_dec = pywt.wavedec(data, pywt.Wavelet(wavelet))
    # Threshold
    wl_thr = map(lambda x: pywt.thresholding.soft(x, threshold), wl_dec)
    # Reconstruct
    return pywt.waverec(wl_thr, wavelet)

#@timer
def wlFlag(data, threshold, wavelet='haar', mode='cpd'):
    """ Wavelet flagging test"""
    (cA, cD) = pywt.dwt(data, wavelet=wavelet, mode=mode)
    
    flags  = np.abs(cD) >= threshold
    cD[flags] = 0
    
    return pywt.idwt(cA, cD, wavelet=wavelet, mode=mode), flags
        

def flagSpectrum(data, thr):
    flags = np.abs(data) >= thr
    return flags

@timer    
def movingAvg(data, window_size):
    """ Apply a moving average to data using convolution"""
    window = np.ones(int(window_size))/float(window_size-1)
    window[window_size/2] = 0 # don't use target spectrum
    
    #window = np.hamming(window_size)
    #print window
    movAvg = convolve1d(data, window, axis=0)
    return movAvg    


@timer
def calibrateSpectrum(freqs, data, tsys):
    """ Calibrate a spectrum using wavelet flagging approach """

    data_w, flags = wlFlag(data_x, 0.005)
    
    tck = interpolate.splrep(freqs[::2][np.invert(flags)], data_w[::2][np.invert(flags)], s=0)
    data_s = interpolate.splev(freqs,tck,der=0)

    data_cal = data/data_s * tsys - tsys

    return data_cal

class sdFits(object):
    def __init__(self, filename):
        fits     = pf.open(filename)
       
        self.filename     = fits.filename()
        self.time    = fits[1].data['TIME']
        self.date    = fits[1].data['DATE-OBS']
        self.data    = fits[1].data['DATA'].astype('float32')
        self.beam    = fits[1].data['BEAM']
        self.flagged = fits[1].data['FLAGGED']
        self.header  = fits[1].header
        self.shape   = np.shape(self.data)
        self.tsys    = fits[1].data['TSYS']

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


freqs    = hi.freqs
data_x   = hi.data[::13,0,0,0,:]
data_y   = hi.data[::13,0,0,1,:]




# Form a RFI free spline for each spectrum
data_avg = movingAvg(data_x, 24*2+1)

for i in range(len(data_avg)):
    data_w, flags = wlFlag(data_avg[i], 0.005)
    tck = interpolate.splrep(freqs[::2][np.invert(flags)], data_w[::2][np.invert(flags)], s=0)
    data_s = interpolate.splev(freqs,tck,der=0)
    data_avg[i] = data_s



tsys     = hi.tsys[::13,0]
tsys  = np.reshape(tsys, [tsys.shape[0],1])
data_cal = data_x / data_avg * tsys - tsys
print data_cal.shape, tsys.shape
print data_cal[32]



plt.subplot(211)
plt.plot(freqs, data_x[33],   color='#333333', label='Raw')
plt.plot(freqs, data_avg[33], color='#cc0000', label='Cal spline')
plt.xlim(1280,1450)
plt.ylim(0,4)
plt.legend()

plt.subplot(212)
plt.plot(freqs, data_x[33]/data_avg[33],   color='#333333', label='Raw / cal')
plt.xlim(1280,1450)
plt.ylim(0.98,1.05)
plt.xlabel("Frequency [MHz]")
plt.legend()

#plt.subplot(212)
#plt.plot(freqs, data_cal[33], color='#333333', label='calibrated')
#plt.xlim(1280,1450)
#plt.ylim(-1,6)

plt.show()

#z = np.abs(freqs-1420)
#ta = z < 1

#plt.plot(freqs[ta], data_cal[33][ta])
#plt.show()
