#!/usr/bin/env python
"""
eager_weaver.py
----------

Functions and helpers to convert HIPSR5 files into SD-FITS

"""

import sys, os, re, time, calendar
from datetime import datetime
import pyfits as pf, numpy as np, tables as tb
from termcolor import cprint

import pylab
from optparse import OptionParser


from lib.printers import LinePrint, Logger
from sdfits import *

__version__  = "v2.0 - Ballistic Bandicoot"
__author__   = "Danny Price"
__email__    = "dprice@cfa.harvard.edu"
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

path = os.getcwd()

def findMatchingTimestamps(h5, sd, gmt_diff=0):
    """ Compare HIPSR and SD-FITS timestamps

    Returns an array of indexes corresponding to the best matches.
    These indexes are for the HDF5 file, and correspond to data rows.
    Throws and error if the t_diff is over 2s (one integration).

    h5: hdf5 file from Hipsr
    sd: sdfits file from mbcorr
    gmt_diff: offset between timestamps in hours. Only required if timestamp is recorded
    in local time, not UTC.
    """

    sd_data = sd[1].data
    hp_ts = h5.root.raw_data.beam_01.col("timestamp")
    hp_dts = np.array([datetime.utcfromtimestamp(ts) for ts in hp_ts])


    utime = sd_data['TIME'][0]
    udate = sd_data['DATE-OBS'][0]

    t_idx = []
    for row in range(len(sd_data['TIME'])):

        utime = sd_data['TIME'][row]
        udate = sd_data['DATE-OBS'][row]

        # From string to datetime obj
        d_d  = datetime.strptime(udate, "%Y-%m-%d")
        # from datetime obj to timestamp
        d_ts = calendar.timegm(d_d.utctimetuple())
        # date + time into timestamp
        dt_ts = d_ts + utime
        # Creating overall timestamp
        dt = datetime.utcfromtimestamp(dt_ts)

        # TODO: Figure out where offset is introduced??!
        t_diffs = hp_ts - dt_ts + gmt_diff * 3600
        idx = np.argmin(np.abs(t_diffs))

        if np.abs(t_diffs[idx]) >= 1.1:
            print "Warning: large t_diff: ",
            print idx, t_diffs[idx]
        t_idx.append(idx)

        if np.abs(t_diffs[idx]) >= 2:
            print "ERROR: Time difference between two files is too large. No match found."
            print "You may have to pass gmt_diff=x, where x is a time offset."
            print "First HIPSR timestamp:  %s"%hp_dts[0]
            utime = sd_data['TIME'][0]
            udate = sd_data['DATE-OBS'][0]
            d_d  = datetime.strptime(udate, "%Y-%m-%d")
            d_ts = calendar.timegm(d_d.utctimetuple())
            dt_ts = d_ts + utime
            dt = datetime.utcfromtimestamp(dt_ts)
            print "First MBCORR timestamp: %s"%dt
            print "Time difference:        %s"%(dt - hp_dts[0])
            exit()

    t_idx = np.array(t_idx)

    return t_idx

def eagerWeaver(sd_file, h5_file, out_file, hp_search_dir=None, sd_search_dir=None, gmt_diff=0):
    """ When trouble comes along you must weave it, weave it good. """

    if hp_search_dir:
        h5_file = os.path.join(hp_search_dir, h5_file)
    if sd_search_dir:
        sd_file = os.path.join(str(sd_search_dir), sd_file)

    print "\nOpening files"
    print "-------------"
    try:
        h5 = tb.openFile(h5_file)
        sd = pf.open(sd_file)
        print "MBCORR: %s"%sd.filename()
        print "HIPSR: %s"%h5.filename


        print "\nGenerating new SD-FITS file"
        print "---------------------------"
        hdulist = generateSDFitsFromMbcorr(sd_file)
        print hdulist

        print "\nMatching timestamps"
        print "-------------------"
        ts_idx = findMatchingTimestamps(h5, sd, gmt_diff)
        print ts_idx

        print "\nFilling in data from HIPSR"
        print "--------------------------"
        pointing = h5.root.pointing.cols
        obs      = h5.root.observation.cols
        sd_data = hdulist[1].data

        print "Rewriting MBCORR common values with HIPSR common values..."
        sd_data["FREQRES"][:]  = np.abs(obs.bandwidth[0])*1e6 / 8192
        sd_data["BANDWID"][:]  = np.abs(obs.bandwidth[0])
        sd_data["CRPIX1"][:]   = 4095
        sd_data["CRVAL1"][:]   = obs.frequency[0] * 1e6
        sd_data["CDELT1"][:]   = np.abs(obs.bandwidth[0])*1e6 / 8192
        sd_data["FLAGGED"][:]  = 0

        scaling = 2**22
        flipped = False
        if obs.bandwidth[0] < 0:
            flipped = True

        # Save file to disk
        skip_on_exist = True
        skip_file = False
        if os.path.exists(out_file):
            if not skip_on_exist:
                print "\nInfo: File %s exists, deleting..."%out_file
                os.remove(out_file)
            else:
                skip_file = True

        if not skip_file:
            print "Filling spectral data...\n"
            for i in range(len(ts_idx)):
                LinePrint("%i of %i"%(i, len(ts_idx)))

                beam_id = "beam_%02d"%sd_data["BEAM"][i]
                h5_row  = ts_idx[i]

                beam = h5.getNode('/raw_data', beam_id)
                xx = beam.cols.xx[h5_row].astype('float32') / scaling
                yy = beam.cols.yy[h5_row].astype('float32') / scaling
                re_xy = beam.cols.re_xy[h5_row].astype('float32') / scaling
                im_xy = beam.cols.im_xy[h5_row].astype('float32') / scaling
                if flipped:
                    xx, yy, re_xy, im_xy = xx[::-1], yy[::-1], re_xy[::-1], im_xy[::-1]
                #data = np.column_stack((xx,yy, re_xy, im_xy))
                data = np.append(xx,yy)
                data = data.reshape([1,1,2,8192])
                sd_data["DATA"][i] = data


            print "\nInfo: Saving to file"
            print out_file
            hdulist.writeto(out_file)
        else:
            print "\nInfo: File %s exists, deleting..."%out_file

        hdulist.close()
        sd.close()
        h5.close()
    except:
        cprint("ERROR:   Could not weave files together.", "red")
        cprint("SD-FITS: %s"%sd_file, "red")
        cprint("HDF5:    %s"%h5_file, "red")
        sd.close()
        h5.close()
        time.sleep(1)
        raise

def findMbcorrFiles(search_dir):
# Regular expression to match hipsr files
    sd_pat = '([0-9A-Za-z-_]+).(sdfits|sdf)'
    filelist = []
    for filename in os.listdir(search_dir):
        match = re.search(sd_pat, filename)
        if match:
            filelist.append(filename)
    return filelist

def findHipsrFiles(search_dir):
# Regular expression to match hipsr files
    hp_pat = '([0-9A-Za-z-_]+).(hdf|h|h5|hdf5)'
    filelist = []
    for filename in os.listdir(search_dir):
        match = re.search(hp_pat, filename)
        if match:
            filelist.append(filename)

    return filelist

def filePairer(sd_filename, search_dir='./'):
    """ Finds closest matching HIPSR file based on SD-FITS filename """

    sd_pat = '(\d+)-(\d+)-(\d+)_(\d\d)(\d\d)-\w+.(sdfits|sdf)'
    match = re.search(sd_pat, sd_filename)

    # Have a first guess that the filenames match exactly
    sd_root = sd_filename.replace(".sdfits", "").replace(".sdf", "")
    hp_root_test = os.path.join(search_dir, sd_root+'.hdf')

    if os.path.exists(hp_root_test):
        return sd_root+'.hdf', 0

    if match:
        # Convert re match to integers, apart from file extension
        (y, m, d, hh, mm) = [int(m) for m in match.groups()[:-1]]
        sd_ts = calendar.timegm((y,m,d,hh,mm,0,0,0,0))
        #print sd_ts

        hp_filelist = findHipsrFiles(search_dir)

        # Version 1 timestamp
        v1_pat = '(P\d+|P\d+s)_(\d+).(hdf|h|h5|hdf5)'
        v2_pat = '(P\d+|P\d+s)_(\d+)-(\d+)-(\d+)_(\d\d)(\d\d)(\d\d).(hdf|h|h5|hdf5)'
        v3_pat = '(\d+)-(\d+)-(\d+)_(\d\d)(\d\d)-([\w_]+).(hdf|h|h5|hdf5)'

        hp_timestamps = []
        for hp_filename in hp_filelist:
            print hp_filename
            v1_match = re.search(v1_pat, hp_filename)
            v2_match = re.search(v2_pat, hp_filename)
            v3_match = re.search(v3_pat, hp_filename)
            if v1_match:
                hp_timestamps.append(int(v1_match.group(2)))
            elif v2_match:
                (y, m, d, hh, mm, ss) = [int(match) for match in v2_match.groups()[1:-1]]
                #print (y,m,d,hh,mm,ss), sd_ts
                hp_timestamps.append(calendar.timegm((y,m,d,hh,mm,ss,0,0,0)))
            elif v3_match:
                (y, m, d, hh, mm) = [int(match) for match in v3_match.groups()[:-2]]
                hp_timestamps.append(calendar.timegm((y,m,d,hh,mm,0,0,0,0)))
            else:
                hp_timestamps.append(0) # to keep array same length as filename array

        hp_timestamps = np.array(hp_timestamps)
        closest = np.argmin(np.abs(hp_timestamps - sd_ts))
        return hp_filelist[closest], np.min(np.abs(hp_timestamps - sd_ts))


if __name__ == '__main__':

    sd_search_dir = '/home/dprice/data/P669-3/mbcorr/sdf'
    output_dir    = '/home/dprice/data/P669-3/hipsr/raw'

    hp_search_dirs = ['/home/dprice/data/P669/hdf/2012-10-29',
                      '/home/dprice/data/P669/hdf/2012-10-30',
                      '/home/dprice/data/P669/hdf/2012-10-31',
                      '/home/dprice/data/P669/hdf/2012-11-01',
                      '/home/dprice/data/P669/hdf/2012-11-02',
                      '/home/dprice/data/P669/hdf/2012-11-03',
                      '/home/dprice/data/P669/hdf/2012-11-04',
                      '/home/dprice/data/P669/hdf/2012-11-05',
                      '/home/dprice/data/P669/hdf/2012-11-06',
                      '/home/dprice/data/P669/hdf/2012-11-07'                  
	               ]

    cprint("HIPSR EAGER WEAVER", 'green')
    cprint("==================\n", 'green')

    time.sleep(0.5)

    print "\nStarting Weave Loop"
    print "-------------------"

    i = 0
    sd_filelist = findMbcorrFiles(sd_search_dir)
    for sd_filename in sd_filelist:
        i += 1
        cprint("\nfile %i of %i (%02d%%)"%(i, len(sd_filelist), float(i)/len(sd_filelist)*100), 'green')
        cprint("-------------------", 'green')

        for hp_search_dir in hp_search_dirs:
            hp_filename, t_diff = filePairer(sd_filename,hp_search_dir)
            if t_diff <= 60:
                break

        if t_diff <= 60:
            print "MBCORR input file:     %s"%sd_filename
            print "Closest matching file: %s"%hp_filename
            print "Time delta: %d\n"%t_diff

            out_filename = os.path.join(output_dir, 'hipsr_'+os.path.basename(sd_filename))
            eagerWeaver(sd_filename, hp_filename, out_filename)
        else:
            print "No matching file found. Skipping..."

