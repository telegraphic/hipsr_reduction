#!/usr/bin/env python
"""
mbcal.py
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


def smoothGauss(list,degree=5):
    """ Gaussian data smoothing
    http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
    """
    window = degree * 2-1
    weight = np.array([1./np.exp( 16. * i * i / window / window) for i in np.arange(-degree + 1, degree)])
    extended = np.r_[[list[0]] * (degree-1), list, [list[-1]] * degree]
    smoothed = np.convolve(weight / weight.sum(), extended, mode='same')

    return smoothed[degree-1:-degree]

def avgDown(col):
    """ Apply average down axis """
    return np.average(col, axis=0).astype('float32')

def flux1934(f):
    """ Return 1934-638 model flux over freq
    frequency in MHz
    """
    log10 = np.log10
    x   = -30.7667 + 26.4908*log10(f) - 7.0977*(log10(f))**2 + 0.605334*(log10(f))**3
    flux =  10**x
    return flux

def flux0407(f):
    """ Return flux of 0407-658
    frequency in MHz
    Note: very rough model (Jan 2014)
    """
    return 6.5 * (f / 2700)**(-1.215)
    
def squash(data, wsize):
    """ Averages together neighbouring bins """
    d = data.reshape((len(data)/wsize, wsize))
    return np.median(d, axis=1)

def mbcal(filename):
    """ Run MBCAL routine to find noise diode calibration """
    print "Opening %s"%filename
    h5 = tb.openFile(filename)

    obs_mode = h5.root.observation[0]["obs_mode"]
    if obs_mode != 'MXCAL':
        print "\nERROR: %s is not an MXCAL file"%filename
        print "ERROR: Please re-run with an MXCAL observation"
        print "ERROR: Observation mode of this file is %s\n"%obs_mode
        exit()

    cal_src = h5.root.pointing[0]["source"]

    if cal_src in ("1934", "1934-638", "PKS B1934-638", "B1934-638"):
        flux_mode = flux1934
    elif cal_src in ("0407", "PKS B0407-658", "0407-658", "B0407-658"):
        flux_model = flux0407
    else:
        print "\nERROR: %s is not a known calibrator"%cal_src
        print "ERROR: Calibration cannot continue!"
        exit()


    # Set up freq axis
    firmware = h5.root.firmware_config[0]["firmware"]
    bw      = h5.root.observation[0]["bandwidth"]
    
    cent_fr = h5.root.observation[0]["frequency"]
    flipped = True if bw < 0 else False
    x_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 8192)
    c_freqs = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, 16)
    c_flux1934 = flux1934(c_freqs)
    f_flux1934 = flux1934(x_freqs)
    if flipped:
        x_freqs, c_freqs, c_flux1934 = x_freqs[::-1], c_freqs[::-1], c_flux1934[::-1]

    # Setup time delta between integrations
    ref_clk   = np.abs(h5.root.observation.cols.bandwidth[0]) * 1e6
    num_chans = h5.root.raw_data.beam_01.cols.xx[0].shape[0]
    acc_len   = h5.root.firmware_config.cols.acc_len[0]
    
    if '400_8192' in firmware:
        # Output is even and odd channels, so factor of 2
        ref_delta = num_chans * acc_len * 2 / ref_clk
    else:
        ref_delta = num_chans * acc_len  / ref_clk

    # Match when beams are pointing at source
    ptime = h5.root.pointing.cols.timestamp[:]
    ids   = h5.root.raw_data.beam_01.cols.id[:]
    ts0   = ptime[0]
    tstamps = [ts0 + (id - ids[0]) * ref_delta for id in ids]
    start_idxs = [np.argmin(np.abs(ptime[ii] - tstamps)) for ii in range(13)]

    
    print "Calibration setup"
    print "-----------------"
    print "Bandwidth:    %s MHz" % bw
    print "Central freq: %s MHz" % cent_fr
    print "Sideband:     %s" % ('lower' if bw < 0 else 'upper')
    print "Number chans: %s" % num_chans
    print "Int. time:    %s s" % ref_delta
    
    
    T_sys_x, T_sys_y = [], []
    for i in range(13):
        beam_id = i+1
        start, stop = start_idxs[i] + 1, start_idxs[i] + 4
                                                                      # Not a typo!
        x_on = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off[start:stop]
        y_on = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off[start:stop]

        try:
            start, stop = start_idxs[i+1] + 1, start_idxs[i+1] + 4
            x_off = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off[start:stop]
            y_off = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off[start:stop]
        except:
            start, stop = start_idxs[i-1] + 1, start_idxs[i-1] + 4
            x_off = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.xx_cal_off[start:stop]
            y_off = h5.getNode("/raw_data/beam_%02i"%beam_id).cols.yy_cal_off[start:stop]
        
        #print beam_id
        #print c_flux1934.shape, x_on.shape, x_off.shape
        #print (avgDown(x_on.astype('float')/x_off) - 1)
        
        Tb_x   = c_flux1934 / (avgDown(x_on.astype('float')/x_off) - 1)
        Tb_y   = c_flux1934 / (avgDown(y_on.astype('float')/y_off) - 1)

        T_sys_x.append(Tb_x)
        T_sys_y.append(Tb_y)

    print "Off-source system temperature (Jy): \n"
    print "    --------------------------"
    print "    | BEAM |  POL A |  POL B |"
    print "    |------------------------|"
    for ii in range(13):
        Tx, Ty = np.array(T_sys_x[ii]), np.array(T_sys_y[ii])
        print "    |  %02i  |  %02.2f |  %02.2f |"%(ii+1, np.average(Tx[4:12]), np.average(Ty[4:12]))
    print "    --------------------------\n"

    temps = np.row_stack((np.array(T_sys_x), np.array(T_sys_y)))

    # Now compute calibration diode temp based upon these data
    cals_x, cals_y = [], []
    ii = 0
    for beam in h5.root.raw_data:

        try:
            start, stop = start_idxs[ii+1] + 1, start_idxs[ii+1] + 4
            c_on_col_x    = beam.col("xx_cal_on").astype('float32')[start:stop]
            c_off_col_x   = beam.col("xx_cal_off").astype('float32')[start:stop]
            c_on_col_y    = beam.col("yy_cal_on").astype('float32')[start:stop]
            c_off_col_y   = beam.col("yy_cal_off").astype('float32')[start:stop]
        except:
            start, stop = start_idxs[ii-1] + 1, start_idxs[ii-1] + 4
            c_on_col_x    = beam.col("xx_cal_on").astype('float32')[start:stop]
            c_off_col_x   = beam.col("xx_cal_off").astype('float32')[start:stop]
            c_on_col_y    = beam.col("yy_cal_on").astype('float32')[start:stop]
            c_off_col_y   = beam.col("yy_cal_off").astype('float32')[start:stop]

        cals_x.append( (avgDown( c_on_col_x / c_off_col_x) - 1) * T_sys_x[ii])
        cals_y.append( (avgDown( c_on_col_y / c_off_col_y) - 1) * T_sys_y[ii])
        ii += 1

    cals_x = np.array(cals_x)
    cals_y = np.array(cals_y)
    cals   = np.row_stack((cals_x, cals_y))

    h5.close()

    print "Noise diode temperatures (Jy): \n"
    print "    ------------------------"
    print "    | BEAM | POL A | POL B |"
    print "    |----------------------|"
    for ii in range(13):
        print "    |  %02i  |  %02.2f |  %02.2f |"%(ii+1, np.average(cals[ii, 4:12]), np.average(cals[ii+13, 4:12]))
    print "    ------------------------\n"

    return cals, temps

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
    print "Saving system temps to %s"%temps_filename
    temps.tofile(temps_filename)