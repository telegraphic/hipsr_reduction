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
import warnings
warnings.filterwarnings("ignore")

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
        temps_filename = os.path.splitext(filename)[0] + '.tsys'
    except:
        print "Usage: mxcal.py [filename]"

    cals, temps = mbcal(filename)

    print "Saving diode calibration to %s"%cal_filename
    cals.tofile(cal_filename)
    #print "Saving system temps to %s"%temps_filename
    #temps.tofile(temps_filename)

    h5 = tb.openFile(filename)
    # Set up freq axis
    bw      = h5.root.observation[0]["bandwidth"]
    cent_fr = h5.root.observation[0]["frequency"]
    flipped = True if bw < 0 else False
    x_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 8192)
    c_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 16)
    c_flux1934 = flux1934(c_freqs)
    f_flux1934 = flux1934(x_freqs)
    if flipped:
        x_freqs, c_freqs, c_flux1934 = x_freqs[::-1], c_freqs[::-1], c_flux1934[::-1]
    h5.close()


    # Plot data

    print "\nPlotting system temperature as function of frequency."
    print "Temperatures should generally lie within dashed lines"
    print "except at the edges of the band (<1200 and >1500 MHz)\n"

    plt.figure(figsize=(12,10))
    plt.subplot(2, 2, 1)
    for i in range(13):
        plt.plot(c_freqs, temps[i], c=colors[i], label='%01i'%(i+1))
    plt.subplot(2, 2, 2)
    for i in range(13):
        plt.plot(c_freqs, temps[13+i], c=colors[i], label='%01i'%(i+1))

    plt.subplot(2, 2, 3)
    for i in range(13):
        plt.plot(c_freqs, cals[i], c=colors[i], label='%01i'%(i+1))
    plt.subplot(2, 2, 4)
    for i in range(13):
        plt.plot(c_freqs, cals[13+i], c=colors[i], label='%01i'%(i+1))

    for i in range(1, 4+1):
        plt.subplot(2, 2, i)
        if i == 1 or i == 3:
            plt.title("Pol A")
        else:
            plt.title("Pol B")

        #plt.xlim(1150,1500)
        if i == 1 or i == 2:
            plt.ylim(0,100)
            plt.ylabel("System temperature [Jy]")
            plt.axhline(y=30, ls='dashed', c='#AAAAAA')
            plt.axhline(y=60, ls='dashed', c='#AAAAAA')
        else:
            plt.ylim(0, 5)
            plt.ylabel("Noise diode temperature [Jy]")
            plt.axhline(y=3.0, ls='dashed', c='#AAAAAA')
            plt.axhline(y=1.5, ls='dashed', c='#AAAAAA')

        plt.xlabel("Frequency [MHz]")
        plt.legend(loc=1)
        plt.minorticks_on()
    plt.tight_layout()
    plt.show()