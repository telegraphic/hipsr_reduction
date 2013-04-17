#!/usr/bin/env python
"""
sdfits.py
----------

Functions and helpers to convert HIPSR5 files into SD-FITS

"""

import sys, os, re, time
from datetime import datetime
import pyfits as pf, numpy as np, tables as tb
import pylab 
from optparse import OptionParser

import lib.config as config
from lib.printers import LinePrint

__version__  = config.__version__
__author__   = config.__author__
__email__    = config.__email__
__license__  = config.__license__
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

path = os.getcwd()

def generateCards(filename):
  """
  Parses a text file and generates a pyfits card list.
  Do NOT feed this a full FITS file, feed it only a human-readable 
  FITS header template. 
  
  A text file is opened, acard is created from each line, then verified. 
  If the line does not pass verification, no card is appended.
  
  Parameters
  ----------
  filename: str
      name of text file header to open and parse
  """
  infile = open(filename)

  header = pf.Header()

  # Loop through each line, converting to a pyfits card
  for line in infile.readlines():
      line = line.rstrip('\n')
      line = line.strip()
      if(line == 'END'):
        break
      else:
        c = pf.Card().fromstring(line)
        c.verify() # This will attempt to fix issuesx[1]
        header.append(c)
        
  return header.cards

def fitsFormatLookup(x):
    """ Helper function to map FITS format codes into numpy format codes
    
    Notes
    -----
    FITS format code         Description                     8-bit bytes
    L                        logical (Boolean)               1
    X                        bit                             *
    B                        Unsigned byte                   1
    I                        16-bit integer                  2
    J                        32-bit integer                  4
    K                        64-bit integer                  4
    A                        character                       1
    E                        single precision floating point 4
    D                        double precision floating point 8
    C                        single precision complex        8
    M                        double precision complex        16
    P                        array descriptor                8    
    """
    
    # This is a python dictionary to lookup mappings.
    return {
           'L' : 'bool_',
           'X' : 'bool_',
           'B' : 'ubyte',
           'I' : 'int16',
           'J' : 'int32',
           'K' : 'int64',
           'A' : 'str_',
           'E' : 'float32',
           'D' : 'float64',
           'C' : 'complex64',
           'M' : 'complex128',
           'P' : 'float32'
    }.get(x, 'float32')

def formatLookup(format_str):
    """ Look up the format of a FITS string """
    pat   = '(\d+)([A-Z])'
    match = re.search(pat, format_str)
    #print match.group()
    
    data_len = int(match.group(1))
    data_fmt = str(match.group(2))
    np_fmt   = fitsFormatLookup(data_fmt)
    np_dtype = '%i%s'%(data_len, np_fmt)
    
    return np_dtype, data_len, np_fmt 

def generateZeros(num_rows, format, dim=None):
    """ Generate blank data to populate binary table 
    
    Used by generateSDFits() to form column definitions.
    
    Parameters
    ----------
    format: str
        FITS format code, e.g. 16A, 2048E
    dim: str
        dimensions for multidimensional data array, e.g. (1024,2,1,1)
        Defaults to None
    """
    
    np_dtype, data_len, np_fmt   = formatLookup(format)
    
    return np.zeros(num_rows, dtype=np_dtype)

def generatePrimaryHDU(hdu_header='header_primaryHDU.txt'):
    """ Generates the Primary HDU
    
    Parameters
    ----------
    hdu_header: string
        Name of the HDU header file to parse to generate the header.
        Defaults to header_primaryHDU.txt
    """
    
    hdu   = pf.PrimaryHDU()
    cards = generateCards(hdu_header)
    
    for card in cards:
        #print card
        if card.keyword == 'COMMENT':
            pass
            hdu.header.add_comment(card.value)
        elif card.keyword == 'HISTORY':
            pass
            hdu.header.add_history(card.value)
        else:
            hdu.header.set(card.keyword, card.value, card.comment)
    
    return hdu

def generateBlankDataHDU(num_rows=1, header_file='header_dataHDU.txt',
                   coldef_file='coldefs_dataHDU.txt'):
    """ Generate a blank data table with N rows.
    
    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_file: str
        Path to the header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions.
        Defaults to 'coldefs_dataHDU.txt'
    
    """
    
    cols = []
    
    # The column definitions are loaded from an external file, which is
    # parsed line-by-line, using regular experssions.
    
    unit_pat   = "unit\s*\=\s*'([\w/%]+)'"
    name_pat   = "name\s*\=\s*'([\w-]+)'"
    dim_pat    = "dim\s*\=\s*'(\([\d,]+\))'"
    format_pat = "format\s*\=\s*'(\w+)'" 

    # Loop through, matching on each line
    cfile = open(coldef_file)
    for line in cfile.readlines():
        unit = name = dim = format = None
        name_match = re.search(name_pat, line)
        if name_match:
            name = name_match.group(1)
             
            format_match = re.search(format_pat, line)
            dim_match    = re.search(dim_pat, line)
            unit_match   = re.search(unit_pat, line)

            if unit_match:   unit = unit_match.group(1)
            if dim_match:    dim  = dim_match.group(1)
                        
            if format_match: 
                fits_fmt = format_match.group(1)
                zarr     = generateZeros(num_rows, fits_fmt, dim)

            
            # Append the column to the column list
            cols.append(pf.Column(name=name, format=fits_fmt, unit=unit, dim=dim, array=zarr))
    
    # Now we have made a list of columns, we can make a new table
    coldefs = pf.ColDefs(cols)
    #print coldefs
    tbhdu   = pf.new_table(coldefs)
    
    # If that all worked, we can populate with the final header values
    cards = generateCards(header_file)
    
    for card in cards:
        if card.key == 'COMMENT':
            pass
            tbhdu.header.add_comment(card.value)
        elif card.key == 'HISTORY':
            pass
            tbhdu.header.add_history(card.value)
        else:
            tbhdu.header.set(card.key, card.value, card.comment)
    
    return tbhdu

def generateDataHDU(input_file, 
                    header_file='lib/header_dataHDU.txt',
                    coldef_file='lib/coldefs_dataHDU.txt'):
    """ Generate a new data table based upon an input file
    
    Parameters
    ----------
    header_file: str
        Path to the header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions.
        Defaults to 'coldefs_dataHDU.txt'
    input_file: str
        String to the input file to grab data from. Defaults to none.
    
    """
    
    sd_in      = pf.open(input_file)
    sd_data    = sd_in[1].data
    num_rows   = sd_data.shape[0]
    
    cols = []
    
    # The column definitions are loaded from an external file, which is
    # parsed line-by-line, using regular experssions.
    
    unit_pat   = "unit\s*\=\s*'([\w/%]+)'"
    name_pat   = "name\s*\=\s*'([\w-]+)'"
    dim_pat    = "dim\s*\=\s*'(\([\d,]+\))'"
    format_pat = "format\s*\=\s*'(\w+)'" 
    
    # Loop through, matching on each line
    cfile = open(coldef_file)
    for line in cfile.readlines():
        unit = name = dim = format = None
        name_match = re.search(name_pat, line)
        if name_match:
            name = name_match.group(1)
             
            format_match = re.search(format_pat, line)
            dim_match    = re.search(dim_pat, line)
            unit_match   = re.search(unit_pat, line)
    
            if unit_match:   
                unit = unit_match.group(1)
            
            
            if dim_match:    
                dim       = dim_match.group(1)
            
            arr_shape = sd_data[name].shape
                    
            if format_match: 
                fits_fmt = format_match.group(1)
                zarr=None

                try:
                    if name == 'DATA' or name == 'FLAGGED':
                        np_dtype, data_len, data_fmt = formatLookup(fits_fmt)
                        print name, " no data"
                    else:
                        # Data array must be flattened (e.g. (2,2) -> 4)
                        np_dtype, data_len, data_fmt = formatLookup(fits_fmt)
                        if data_len > 1 and data_fmt != 'str_':
                            z_shape = (sd_data[name].shape[0], data_len)
                        else:
                             z_shape = sd_data[name].shape
                        #print name, z_shape, sd_data[name].shape
                        zarr     = sd_data[name].reshape(z_shape)
                        
                except:
                    print "Error with %s"%name
            
                # Append the column to the column list
                cols.append(pf.Column(name=name, format=fits_fmt, unit=unit, dim=dim, array=zarr))
    
    # Now we have made a list of columns, we can make a new table
    #print cols
    coldefs = pf.ColDefs(cols)
    #print coldefs
    tbhdu   = pf.new_table(coldefs)
    
    # If that all worked, we can populate with the final header values
    cards = generateCards(header_file)
    
    for card in cards:
        if card.keyword == 'COMMENT':
            pass
            tbhdu.header.add_comment(card.value)
        elif card.keyword == 'HISTORY':
            pass
            tbhdu.header.add_history(card.value)
        else:
            tbhdu.header.set(card.keyword, card.value, card.comment)
    
    return tbhdu


def timestamp2dt(timestamp):
    """ Convert timestamp to date and time for SD-FITS """
    
    dt = datetime.utcfromtimestamp(timestamp)
    
    date = dt.strftime("%Y-%m-%d")
    # TODO: Check this is correct
    time = dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond * 1e-6
    return (date, time)
    
def findMatchingTimestamps(h5, sd):
    """ Compare HIPSR and SD-FITS timestamps

    Returns an array of indexes corresponding to the best matches.
    These indexes are for the HDF5 file, and correspond to data rows.
    Throws and error if the t_diff is over 2s (one integration).

    """

    sd_data = sd[1].data
    hp_ts = h5.root.raw_data.beam_01.col("timestamp")
    hp_dts = np.array([datetime.utcfromtimestamp(ts) for ts in hp_ts])
    
    t_idx = []
    for row in range(len(sd_data['TIME'])):
            
        utime = sd_data['TIME'][row]
        udate = sd_data['DATE-OBS'][row]
        
        # From string to datetime obj
        d_d  = datetime.strptime(udate, "%Y-%m-%d")
        # from datetime obj to timestamp
        d_ts = time.mktime(d_d.utctimetuple())
        # date + time into timestamp
        dt_ts = d_ts + utime
        # Creating overall timestamp
        dt = datetime.utcfromtimestamp(dt_ts)

        # TODO: Figure out why an 8hr offset is required??!
        t_diffs = hp_ts - dt_ts - 60*60*8
        idx = np.argmin(np.abs(t_diffs))

        if np.abs(t_diffs[idx]) >= 1.1:
            print "Warning: large t_diff: ",
            print idx, t_diffs[idx]
        t_idx.append(idx)

        if np.abs(t_diffs[idx]) >= 2:
            print "ERROR: Time difference between two files is too large. No match found."
            raise

    t_idx = np.array(t_idx)

    return t_idx

def generateBlankSDFits(num_rows,
                        header_primary='lib/header_primaryHDU.txt',
                        header_tbl='lib/header_dataHDU.txt',
                        coldef_file='lib/coldefs_dataHDU.txt'):
    """ Generate a blank SD-FITS file

    This function returns a blank SD-FITS file, with num_rows in the binary table.
    It generates all the required columns, then fills them with blank data (zeros).

    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_primary: str
            Path to the primaryHDU header file. Defaults to 'header_primaryHDU.txt'
    header_data: str
        Path to the binary table header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions for the binary table.
        Defaults to 'coldefs_dataHDU.txt'

    """

    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateBlankDataHDU(num_rows, header_tbl, coldef_file)
    hdulist = pf.HDUList([prhdu, tbhdu])

    return hdulist

def generateSDFitsFromMbcorr(input_file,
                             header_primary='lib/header_primaryHDU.txt',
                             header_tbl='lib/header_dataHDU.txt',
                             coldef_file='lib/coldefs_dataHDU.txt'):
    """ Generate a blank SD-FITS file

    This function returns a blank SD-FITS file, with num_rows in the binary table.
    It generates all the required columns, then fills them with blank data (zeros).

    Parameters
    ----------
    num_rows: int
        The number of rows in the binary table.
    header_primary: str
            Path to the primaryHDU header file. Defaults to 'header_primaryHDU.txt'
    header_data: str
        Path to the binary table header file. Defaults to 'header_dataHDU.txt'
    coldef_file: str
        Path to the file containing column definitions for the binary table.
        Defaults to 'coldefs_dataHDU.txt'

    """

    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateDataHDU(input_file, header_tbl, coldef_file)
    hdulist = pf.HDUList([prhdu, tbhdu])

    return hdulist

def generateSDFitsFromHipsr(filename_in, path_in, filename_out, path_out):
    """ Generate an SD-FITS file from a hipsr5 file """
    
    # Open h5 file
    h5file = os.path.join(path_in, filename_in)
    h5 = tb.openFile(h5file)
    
    num_acc  = h5.root.raw_data.beam_01.shape[0] 
    num_rows = num_acc * 13
    
    if num_acc == 0:
        print "No data in %s. Skipping."%h5file
        return -1
    
    print "Input file: %s"%h5.filename
    print "No accumulations: %s, no rows: %s"%(num_acc, num_rows)
    
    # We now need to generate a blank SD-FITS file, with the same number of rows
    print "\nGenerating blank SD-FITS file with %i rows..."%num_rows
    hdulist = generateBlankSDFits(num_rows)
    print hdulist.info()
    
    # Next, we copy over observation data    
    print "Filling new SD-FITS with HIPSR data..."
    
    pointing = h5.root.pointing.cols
    obs      = h5.root.observation.cols
    sdtab    = hdulist[1].data
    
    # Fill in common values
    print "Filling in common values... ",
    sdtab["SCAN"][:]     = 1
    sdtab["EXPOSURE"][:] = obs.acc_len[0]
    sdtab["OBJECT"][:]   = pointing.source[0]
    sdtab["OBJ-RA"][:]   = pointing.ra[0]
    sdtab["OBJ-DEC"][:]  = pointing.dec[0]
    sdtab["RESTFRQ"][:]  = obs.frequency[0]    
    sdtab["FREQRES"][:]  = np.abs(obs.bandwidth[0])*1e6 / 8192
    sdtab["BANDWID"][:]  = np.abs(obs.bandwidth[0])
    sdtab["CRPIX1"][:]   = 4095
    sdtab["CRVAL1"][:]   = obs.frequency[0] * 1e6
    sdtab["CDELT1"][:]   = np.abs(obs.bandwidth[0])*1e6 / 8192
    sdtab["FLAGGED"][:]  = 0
    sdtab["SCANRATE"][:] = obs.scan_rate[0]


    # No TCS info - common, TODO ASAP
    sdtab["OBSMODE"][:]  = 'SC' 
    sdtab["IF"][:]       = 1
    print "OK."
    
    row_sd   = 0
    cycle_id = 0
    scaling = 2**22 # Divide through to change 32-bit to 
    
    flipped = False
    if obs.bandwidth[0] < 0:
        flipped = True
    
    print "Filling in unique values... "
    scan_pointing_len = h5.root.scan_pointing.shape[0]
    
    for row_h5 in range(num_acc):
        cycle_id += 1 # Starts at 1 in SD-FITS file
        for beam in h5.root.raw_data:
            LinePrint("%i of %i"%(row_sd, num_rows))
            
            if cycle_id <= scan_pointing_len:
                raj_id = "mb%s_raj"%beam.name.lstrip('beam_')
                dcj_id = "mb%s_dcj"%beam.name.lstrip('beam_')
                
                sdtab["CYCLE"][row_sd]   = cycle_id
                beam_id = int(beam.name.lstrip('beam_'))
                
                # Fix beam mapping (remove after fixing mapping)
                if beam_id == 11:   beam_id = 13
                elif beam_id == 13: beam_id = 11
                
                sdtab["BEAM"][row_sd]     = beam_id
                
                sdtab["CRVAL3"][row_sd]   = h5.root.scan_pointing.col(raj_id)[cycle_id-1]
                sdtab["CRVAL4"][row_sd]   = h5.root.scan_pointing.col(dcj_id)[cycle_id-1]
                sdtab["AZIMUTH"][row_sd]  = h5.root.scan_pointing.col("azimuth")[cycle_id-1]
                sdtab["ELEVATIO"][row_sd] = h5.root.scan_pointing.col("elevation")[cycle_id-1]
                
                try:
                    timestamp  = beam.cols.timestamp[row_h5]
                    date_obs, time = timestamp2dt(timestamp)
                    sdtab["DATE-OBS"][row_sd] = date_obs
                    sdtab["TIME"][row_sd]     = time
                    
                    xx = beam.cols.xx[row_h5].astype('float32') / scaling
                    yy = beam.cols.yy[row_h5].astype('float32') / scaling
                    re_xy = beam.cols.re_xy[row_h5].astype('float32') / scaling
                    im_xy = beam.cols.im_xy[row_h5].astype('float32') / scaling
                
                    if flipped:
                        xx, yy, re_xy, im_xy = xx[::-1], yy[::-1], re_xy[::-1], im_xy[::-1]
                    
                    data = np.column_stack((xx,yy, re_xy, im_xy))
                    data = data.reshape([1,1,4,8192])    
                    sdtab["DATA"][row_sd] = data
                    
                except:
                    print "\nWARNING: missing row in %s"%beam.name
                    print "Current index: %i"%row_h5
                    print "Row length: %i"%beam.shape[0]
                

                row_sd += 1
            else:
                print "WARNING: scan_pointing table is not complete."
                print "%s table length: %i"%(beam.name, beam.shape[0])
                print "scan_pointing table length: %i"%scan_pointing_len
    
    print "\nWriting to file...",
    hdulist.writeto(os.path.join(path_out, filename_out))
    h5.close()
    print "OK"


if __name__ == '__main__':

    # Setup input and output files
    sd_file  = '/media/BACKUP_DATA/HIPSR/P669/mbcorr/sdf/2012-11-03_1226-P669_g02_1285_p669_dec01a.sdfits'
    h5_file  = '/media/BACKUP_DATA/HIPSR/P669/2012-11-03/2012-11-03_1226-P669_g02_1285_p669_dec01a.hdf'
    out_file = '/home/dprice/hipsr/scratch/2012-11-03_1226-P669_g02_1285_p669_dec01a.sdfits'

    print "\nOpening files"
    print "-------------"
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
    ts_idx = findMatchingTimestamps(h5, sd)
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
        data = np.column_stack((xx,yy))
        data = data.reshape([1,1,2,8192])
        sd_data["DATA"][i] = data


    # Save file to disk
    if os.path.exists(out_file):
        print "\nInfo: File exists, deleting..."
        os.remove(out_file)

    print "\nInfo: Saving to file"
    hdulist.writeto(out_file)
    hdulist.close()
