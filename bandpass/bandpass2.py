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

@timer
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
    window = np.ones(int(window_size))/float(window_size)
    movAvg = convolve1d(d, window, axis=0)
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

data_x    = hi.data[0,0,0,0,:].astype('float32')
freqs     = hi.freqs
data_w, flags = wlFlag(data_x, 0.005)
tsys = hi.tsys[0,0]


tck = interpolate.splrep(freqs[::2][np.invert(flags)], data_w[::2][np.invert(flags)], s=0.001)
data_s = interpolate.splev(freqs,tck,der=0)
f = freqs[::2][np.invert(flags)]
d = data_x[::2][np.invert(flags)]

data_cal = data_x/data_s * tsys - tsys
std = np.std(data_cal[::2][np.invert(flags)])
print "St.Dev of calibrated data: %s"%std

data_cal = calibrateSpectrum(freqs, data_x, tsys)
flags    = np.invert(flagSpectrum(data_cal, 5*std))

std = np.std(data_cal[flags][0:len(data_cal[flags])/2])
print "St.Dev of low freq data: %s"%std
std = np.std(data_cal[flags][len(data_cal[flags])/2:])
print "St.Dev of high freq data: %s"%std
std = np.std(data_cal[flags][len(data_cal[flags])/4:3*len(data_cal[flags])/4])
print "St.Dev of mid freq data: %s"%std

data_all = hi.data[:,0,0,0,:].astype('float32')
print data_all.shape
data_avg = movingAvg(data_x, 32)


plt.subplot(211)
plt.plot(freqs, data_x, color='#333333', label='Original data')
plt.plot(f, d, color='#cc0000', label='Post flagging')
plt.plot(freqs, data_s, color='#006600', label='Splined')
plt.xlim(1150,1500)
plt.ylim(-0.09,6)
plt.legend()


plt.subplot(212)
plt.plot(freqs, data_cal, color='#333333', label='Calibrated data')
plt.plot(freqs[flags], data_cal[flags], color='#CC0000', label='Flagged data')
plt.legend()
plt.xlim(1150,1500)
plt.ylim(-10*std, 100*std)
plt.show()

