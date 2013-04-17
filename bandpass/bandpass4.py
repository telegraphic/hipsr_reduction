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
    window[window_size/2] = 0
    #window = np.hamming(window_size)
    print window
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
        self.data    = fits[1].data['DATA']
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
data_x   = hi.data[:,0,0,0,:].astype('float32')
data_y   = hi.data[:,0,0,1,:].astype('float32')
tsys     = hi.tsys[0,0]


data_avg = movingAvg(data_x, 49)

# Form a RFI free spline for each spectrum
for i in range(len(data_avg)):
    data_w, flags = wlFlag(data_avg[i], 0.005)
    tck = interpolate.splrep(freqs[::2][np.invert(flags)], data_w[::2][np.invert(flags)], s=0)
    data_s = interpolate.splev(freqs,tck,der=0)
    data_avg[i] = data_s

tsys     = hi.tsys[:,0]
tsys  = np.reshape(tsys, [tsys.shape[0],1])
data_cal = data_x / data_avg * tsys - tsys
print data_cal.shape, tsys.shape
print data_cal[32]



'''plt.subplot(211)
plt.plot(freqs, data_x[32],   color='#333333', label='Raw')
plt.plot(freqs, data_avg[32], color='#cc0000', label='Spline')
plt.xlim(1200,1450)
plt.ylim(0,10)
plt.legend()'''

d = data_cal[32]
wl = 'db2'
mode = 'cpd'

(cA0, cD0)   = pywt.dwt(d, wl, mode)
print "level 0", len(cA0), len(cD0)

(cA1, cD1) = pywt.dwt(cA0, wl, mode)
print "level 1", len(cA1), len(cD1)

(cA2, cD2) = pywt.dwt(cA1, wl, mode)
print "level 2", len(cA2), len(cD2)

(cA3, cD3) = pywt.dwt(cA2, wl, mode)
print "level 3", len(cA3), len(cD3)

(cA4, cD4) = pywt.dwt(cA3, wl, mode)
print "level 4", len(cA4), len(cD4)

(cA5, cD5) = pywt.dwt(cA4, wl, mode)
print "level 5", len(cA5), len(cD5)

(cA6, cD6) = pywt.dwt(cA5, wl, mode)
print "level 6", len(cA6), len(cD6)

plt.subplot(421)
plt.plot(cA0)
plt.ylim(-20,50)

plt.subplot(422)
plt.plot(cD0)
plt.ylim(-5,5)

plt.subplot(423)
plt.plot(cA1)
plt.ylim(-20,50)

plt.subplot(424)
plt.plot(cD1)
plt.ylim(-5,5)

plt.subplot(425)
plt.plot(cA2)
plt.ylim(-20,50)

plt.subplot(426)
plt.plot(cD2)
plt.ylim(-5,5)

plt.subplot(427)
plt.plot(cA3)
plt.ylim(-20,50)

plt.subplot(428)
plt.plot(cD3)
plt.ylim(-5,5)

plt.show()

plt.plot(freqs, d)
#cA6 = np.zeros_like(cA6)

thr = 1
cA3 = np.zeros_like(cA3)
cD2[cD2 < thr] = 0
cD1[cD1 < thr] = 0

#print "level 6", len(cA5), len(cD5)
#cA5 = pywt.idwt(cA6, cD6, wl)
#print "level 5", len(cA5), len(cD5)
#cA4 = pywt.idwt(cA5, cD5, wl)
#print "level 4", len(cA4), len(cD4)
#cA3 = pywt.idwt(cA4, cD4, wl, mode)
print "level 3", len(cA3), len(cD3)
cA2 = pywt.idwt(cA3, cD3, wl, mode)
print "level 2", len(cA2), len(cD2)
cA1 = pywt.idwt(cA2[:len(cD2)], cD2, wl, mode)
print "level 1", len(cA1), len(cD1)
cA0 = pywt.idwt(cA1[:len(cD1)], cD1, wl, mode)
print "level 0", len(cA0), len(cD0)
dw  = pywt.idwt(cA0[:len(cD0)], cD0, wl, mode)
plt.plot(freqs, dw)
plt.ylim(-10, 25)
plt.show()

