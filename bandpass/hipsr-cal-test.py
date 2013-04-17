#!/usr/bin/env python

"""
hipsr-calibrate.py
==================

Bandpass calibration and RFI flagging for HIPSR data.

This script implements Barnes et. al. (2001) calibration:

S = I / B * T_bp - T_tar

S     = calibrated spectrum
I     = uncalibrated spectrum
B     = bandpass estimate
T_bp  = Bandpass system temperature (off target)
T_tar = System temp while on target

"""

import sys, os, re, time
import matplotlib, pylab as plt

import numpy as np
from scipy import interpolate
from scipy.ndimage.filters import convolve1d

import pyfits as pf
import pywt


def timer(func):
    """ Decorator for timing function calls"""
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper
    
class LinePrint():
    """
    Print things to stdout on one line dynamically
    """
    def __init__(self,data):
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()    

def wlFlag(data, threshold, wavelet='haar', mode='cpd'):
    """ Wavelet flagging subroutine """
    (cA, cD) = pywt.dwt(data, wavelet=wavelet, mode=mode)
    flags  = np.abs(cD) >= threshold
    cD[flags] = 0
    return pywt.idwt(cA, cD, wavelet=wavelet, mode=mode), flags


def movingAvg(data, window_size):
    """ Apply a moving average to data using convolution"""
    
    #window = np.ones(int(window_size))/float(window_size-1)
    #window[window_size/2] = 0 # don't use target spectrum
    
    # Using Hanning window to weight data closer higher
    window   = np.hanning(window_size)
    window   = window / window.sum()
    window   = np.insert(window, len(window)/2, 0) 
    
    movAvg = convolve1d(data, window, axis=0)
    return movAvg    

def spliner(data, freqs):
    """Helper fn for toSpine2"""
    data_w, flagged = wlFlag(data, 0.005)
    unflagged = np.invert(flagged)
    unflagged = np.ones_like(flagged)
    tck = interpolate.splrep(freqs[::2][unflagged], data[::2][unflagged], s=0)
    data_s = interpolate.splev(freqs, tck, der=0)
    return data    

@timer
def toSpline2(freqs, data):
    """ Thing that I'm going to vectorize"""
    return np.apply_along_axis(spliner, 1, data, freqs)

@timer
def toSpline(freqs, data):
    """ Creates calibration splines from data """
    for i in range(len(data)): 
        data_w, flagged = wlFlag(data[i], 0.005)
        unflagged = np.invert(flagged)
        tck = interpolate.splrep(freqs[::2][unflagged], data_w[::2][unflagged], s=1)
        data_s = interpolate.splev(freqs, tck, der=0)
        data[i] = data_s
    return data


def applyCal(raw_data, cal_data, t_bp, t_tar):
    """ Apply calibration to data 
    
    Using Barnes et. al. (2001) eqn. 3:
    S = I / B * T_bp - T_tar

    S     = calibrated spectrum
    I     = uncalibrated spectrum
    B     = bandpass estimate
    T_bp  = Bandpass system temperature (off target)
    T_tar = System temp while on target
    
    """
    return (raw_data / cal_data * t_bp) - t_tar


def flagSpectrum(spectra, sigma):
    """ Flag spectral features about a given RMS sigma """
    l = len(spectra[0])
    
    spectra_midband = spectra[:, l/2-l/8:l/2+l/8]
    med = np.median(np.abs(spectra_midband), axis=1)
    
    stds, avgs = np.zeros((len(med),1)), np.zeros((len(med),1))
    for row in range(len(spectra_midband)):
        sm, mrow = spectra_midband[row], med[row]
        spectra_bel_med = sm[sm < sigma * mrow]
        
        stds[row,0] = np.std(spectra_bel_med)
        avgs[row,0] = np.average(spectra_bel_med)
    
    flagged = spectra > 10
    for row in range(len(spectra)):
        flagged[row] = np.abs(spectra[row] - avgs[row]) >= stds[row] * sigma
    
    return flagged

#@timer
def removeResidualFlux(freqs, spectra, flags, poly_order):
    """ Remove residual flux by fitting an nth order polynomial"""
    
    for row in range(len(spectra)):
        unflg = flags[row] == 0
        fit  = np.polyfit(freqs[unflg], spectra[row][unflg], poly_order)
        poly = np.poly1d(fit)
        spectra[row] = spectra[row] - poly(freqs)
    
    return spectra

class sdFits(object):
    """ SD-FITS helper class """
    def __init__(self, filename):
        fits     = pf.open(filename)
        
        self.fits    = fits
        self.filename= fits.filename()
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

@timer
def calibrate(filename_in, filename_out):

    hi = sdFits(filename_in)
    #print "HIPSR:  %s"%hi.filename

    cycles = 12*2+1
    sigma  = 3
    data_glob    = np.copy(hi.data[:,:,:,:,:])
    flags_glob   = np.copy(hi.flagged)
    
    #print "Calibrating spectra..."
    # Create processes    
    for i in range(0,13):
        LinePrint("Processing beam %i of %i"%(i+1, 13))
        freqs    = hi.freqs
        data_a   = hi.data[i::13,0,0,0,:]
        data_b   = hi.data[i::13,0,0,1,:]
        tsys_a   = hi.tsys[i::13,0]
        tsys_b   = hi.tsys[i::13,1]
        
        # Reshape t_sys
        t_sys_a = np.reshape(tsys_a, [tsys_a.shape[0],1])
        t_sys_b = np.reshape(tsys_b, [tsys_b.shape[0],1])

        # Compute the rolling average of each spectrum
        data_avg_a = movingAvg(data_a, cycles)
        data_avg_b = movingAvg(data_b, cycles)

        # Convert average spectra to RFI-free splines
        data_s_a   = toSpline2(freqs, data_avg_a)
        data_s_b   = toSpline2(freqs, data_avg_b)
        
        # Compute T_bp using rolling average
        t_bp_a    = movingAvg(t_sys_a, cycles)
        t_bp_b    = movingAvg(t_sys_b, cycles) 
        
        # Apply calibration
        data_cal_a = applyCal(data_a, data_s_a, t_bp_a, t_sys_a)
        data_cal_b = applyCal(data_b, data_s_b, t_bp_b, t_sys_b)
        
        # Flag spectrum
        flags_a = flagSpectrum(data_cal_a, sigma)
        flags_b = flagSpectrum(data_cal_b, sigma)
        
        flags_a[:, freqs > 1480] = 1
        flags_b[:, freqs > 1480] = 1

        flags_a[:, freqs < 1284] = 1
        flags_b[:, freqs < 1284] = 1       

        # Apply cal correction to set flux to 0Jy avg
        data_cal_a = removeResidualFlux(freqs, data_cal_a, flags_a, 1)
        data_cal_b = removeResidualFlux(freqs, data_cal_b, flags_a, 1)

        # Unflag rest HI line       
        flags_a[:, np.abs(freqs - 1420.4) < 0.5] = 0
        flags_b[:, np.abs(freqs - 1420.4) < 0.5] = 0
        
        flags_a[:, np.abs(freqs - 1450.5) < 1.0] = 1
        flags_b[:, np.abs(freqs - 1450.5) < 1.0] = 1

        # Write to new data array
        data_glob[i::13,0,0,0,:]  = data_cal_a
        data_glob[i::13,0,0,1,:]  = data_cal_b
        flags_glob[i::13,0,0,0,:] = flags_a
        flags_glob[i::13,0,0,1,:] = flags_b
       
    # Write to new file
    print "\nWriting to new file..."
    cal_fits = hi.fits
    d_len    = len(cal_fits[1].data["DATA"])
    for row in range(d_len):
        LinePrint("writing %i or %i"%(row, d_len))
        cal_fits[1].data["DATA"][row]    = data_glob[row]
        cal_fits[1].data["FLAGGED"][row] = flags_glob[row]
        
    print "\nSaving to %s..."%filename_out
    #if os.path.exists(filename_out): os.remove(filename_out)
    cal_fits.writeto(filename_out)
    
    
if __name__ == '__main__':

    filepath_in  = '/home/share/data/P669/hipsr/raw/hipsr_2012-10-29_1113-P669_g02_1285_p669_dec02a.sdfits'
    filepath_out = './test2.sdfits'
    calibrate(filepath_in, filepath_out)
  
    
    
