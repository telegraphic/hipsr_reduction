#!/usr/bin/env python
"""
hipsr-mbcal.py
==============

Computes the noise diode temperature in Jy from observation of a calibration source.

This script returns the temperature of the noise diode for each beam and polarization,
as a function of frequency (16x25 MHz channels over 400 MHz).

Notes
-----

Calibration is done using a Y-factor measurement, based on an assumed model of the calibration
source. We take four measurements per beam:

    P_src:  Power when on-source, with noise diode off
    P_bg :  Power when off-source, with noise diode off
    P_Don:  Power when off-source, with noise diode on
    P_Doff: Power when off-source, with noise diode off == P_bg

And convert these into temperatures (in Jy, not K!)

    T_sys: System temperature of telescope when off-source, diode off
    T_D:   Equivalent noise temperature of diode

The relations being applied are:

    T_sys  = T_src / (P_src / P_bg - 1)
    T_D    = (P_Don / P_Doff - 1) * T_sys

Once T_D has been measured, we use this to compute the Tsys for the telescope when we observe
sources of interest:

    T_sys' = T_D / (P_don / P_doff - 1)

The output of this script is a (26 x 16) numpy array of T_D as a function of frequency, which
can be used for calibrating future data.

"""
import sys, os
import numpy as np
import tables as tb

import pylab as plt

params = {'legend.fontsize': 10,
              'legend.linewidth': 1.5}
plt.rcParams.update(params)
from datetime import datetime
from lib.hipsrx import Hipsr6

from lib.mbcal import *

colors = [
    '#cd4a4a', '#ff6e4a', '#9f8170', '#ffcf48', '#bab86c', '#c5e384', '#1dacd6',
    '#71bc78', '#9aceeb', '#1a4876', '#9d81ba', '#cdc5c2', '#fc89ac'
    ]

if __name__ == '__main__':
    
    # Open mxcal file
    try:
        filename = sys.argv[1]
        cal_filename = os.path.splitext(filename)[0] + '.cal'
    except:
        print "Usage: mxcal.py [filename]"
    
    cals = mbcal(filename)

    print "Saving diode calibration to %s"%cal_filename
    cals.tofile(cal_filename)
    
    #print cals.shape
    #print "Freqs:      %s"%['%2.2f'%v for v in c_freqs[5:9]]
    #print "Average:    %s"%['%2.2f'%v for v in np.average(cals_y[:,5:9], axis=0)]
    #print "Median:     %s"%['%2.2f'%v for v in np.median(cals_y[:,5:9], axis=0)]
    #print "Stdev:      %s"%['%2.2f'%v for v in np.std(cals_y[:,5:9], axis=0)]
    #print "Fractional: %s"%['%2.2f'%v for v in (np.std(cals_y[:,5:9], axis=0) / np.average(cals_y[:,4:8], axis=0) * 100)]
     
    #plt.subplot(212)
    #for i in range(13):
    #    plt.plot(c_freqs, cals[i], c=colors[i], label='%01i'%(i+1))
    ##plt.xlim(1150,1500)
    #plt.ylim(0,5)
    #plt.xlabel("Frequency [MHz]")
    #plt.ylabel("Noise diode temp [Jy]")
    #plt.legend(loc=2, ncol=7)
    #plt.tight_layout()
    #plt.show()

