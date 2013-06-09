#!/usr/bin/env python
"""
cal_diode_temp.py
-----------

First attempt to calibrate using the NAR data

"""
import sys, os
import numpy as np
import tables as tb
import pyfits as pf
import pylab as plt
from datetime import datetime
from lib.hipsr6 import Hipsr6

def avgDown(col):
    """ Apply average down axis """
    return np.median(col, axis=0).astype('float32')

def flux1934(f):
    """ Return 1934-638 model flux over freq
    frequency in MHz
    """
    log10 = np.log10
    x   = -30.7667 + 26.4908*log10(f) - 7.0977*(log10(f))**2 + 0.605334*(log10(f))**3
    flux =  10**x
    return flux
    
def squash(data, wsize):
    """ Averages together neighbouring bins """
    
    d = data.reshape((len(data)/wsize, wsize))
    return np.median(d, axis=1)

if __name__ == '__main__':
    
    colors = [
        '#cd4a4a', '#ff6e4a', '#9f8170', '#ffcf48', '#bab86c', '#c5e384', '#1dacd6',
        '#71bc78', '#9aceeb', '#1a4876', '#9d81ba', '#cdc5c2', '#fc89ac'
        ]

    file_path = 'obs_1934'
    file_list = os.listdir(file_path)
    file_list.sort()
    
    h = []
    for f in file_list:
        b = tb.openFile(os.path.join(file_path, f))
        print "Beam on source: %i"%b.root.observation[0]["ref_beam"]
        h.append(b)
    
    # Set up freq axis
    bw      = h[0].root.observation[0]["bandwidth"]
    cent_fr = h[0].root.observation[0]["frequency"]
    flipped = True if bw < 0 else False
    
    x_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 8192)
    c_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 16)
    c_flux1934 = flux1934(c_freqs)
    
    if flipped:
        print "Flipped"
        x_freqs, c_freqs = x_freqs[::-1], c_freqs[::-1]
    
    # Compute T_sys from ON / OFF source
    T_sys_x, T_sys_y = [], []
    for i in range(len(h)):
        beam_id = i+1
                                                                      # Not a typo!
        x_on = avgDown(h[i].getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off)
        y_on = avgDown(h[i].getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off)
        try:
            x_off = avgDown(h[i+1].getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off)
            y_off = avgDown(h[i+1].getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off)
        except:
            x_off = avgDown(h[i-1].getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off)
            y_off = avgDown(h[i-1].getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off)
        
        Tb_x   = c_flux1934 / (x_on/x_off -1)
        Tb_y = c_flux1934 / (y_on/y_off -1)  
        
        T_sys_x.append(Tb_x)
        T_sys_y.append(Tb_y)
        plt.plot(c_freqs, Tb_y, label=beam_id, c=colors[i])
    
    plt.ylim(0,100)
    plt.xlim(1150,1500)
    plt.legend(loc=3)
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("System temp [Jy]")
    #plt.savefig('plots_parkes/tsys_13_beams.pdf')
    plt.show()
    
    
    cals_x, cals_y = [], []
    for beam in h[0].root.raw_data:        
        c_on_col_x    = beam.col("xx_cal_on").astype('float32')
        c_off_col_x   = beam.col("xx_cal_off").astype('float32')
        c_on_col_y    = beam.col("yy_cal_on").astype('float32')
        c_off_col_y   = beam.col("yy_cal_off").astype('float32')
        
        cals_x.append( (avgDown( c_on_col_x / c_off_col_x) - 1) * T_sys_x[i])
        cals_y.append( (avgDown( c_on_col_y / c_off_col_y) - 1) * T_sys_y[i])

    # Override beam 01 as it's on source! 
    beam = h[1].root.raw_data.beam_01     
    c_on_col_x    = beam.col("xx_cal_on").astype('float32')
    c_off_col_x   = beam.col("xx_cal_off").astype('float32')
    c_on_col_y    = beam.col("yy_cal_on").astype('float32')
    c_off_col_y   = beam.col("yy_cal_off").astype('float32')
        
    cals_x[0] = (avgDown( c_on_col_x / c_off_col_x) - 1) * T_sys_x[0] 
    cals_y[0] = (avgDown( c_on_col_y / c_off_col_y) - 1) * T_sys_y[0] 
    
    cals_x = np.array(cals_x)
    cals_y = np.array(cals_y)
    
    #np.savetxt('diode_temps_Apr_2013.csv', cals)
    print "Saving cals"
    cals_x.tofile('diode_jy_x.cal')
    cals_y.tofile('diode_jy_y.cal')
    
    print "Freqs:      %s"%['%2.2f'%v for v in c_freqs[5:9]]
    print "Average:    %s"%['%2.2f'%v for v in np.average(cals_y[:,5:9], axis=0)]
    print "Median:     %s"%['%2.2f'%v for v in np.median(cals_y[:,5:9], axis=0)]
    print "Stdev:      %s"%['%2.2f'%v for v in np.std(cals_y[:,5:9], axis=0)]
    print "Fractional: %s"%['%2.2f'%v for v in (np.std(cals_y[:,5:9], axis=0) / np.average(cals_y[:,4:8], axis=0) * 100)]
     
    # Plot the noise diode cal  
    for i in range(len(cals_y)):
        plt.plot(c_freqs, cals_y[i], c=colors[i], label='%01i'%(i+1))
    plt.xlim(1100,1500)
    plt.ylim(0,5)
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Noise diode temp [Jy]")
    plt.legend(loc=3)
    #plt.savefig('nd_temp_13_beams.pdf')
    plt.show()
    
    for hdf in h:
        hdf.close()
