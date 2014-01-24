#!/usr/bin/env python
"""
sdfits.py
----------

Functions and helpers to convert HIPSR5 files into SD-FITS

"""

import sys, os, re, time
from datetime import datetime
import numpy as np, tables as tb

from lib.printers import LinePrint
from lib.hipsrx import Hipsr6
from lib.mbcal import mbcal

try:
    import pyfits as pf
except ImportError:
    try:
        from astropy.io import fits as pf
        print "Using Astropy for FITS I/O"
    except:
        print "Error: cannot load PyFITS or AstroPY I/O. Please check your install."
        exit()

__version__  = "v2.0 - Ballistic Bandicoot"
__author__   = "Danny Price"
__email__    = "dprice@cfa.harvard.edu"
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

path = os.path.abspath(__file__).replace('sdfits.pyc', '').replace('sdfits.py', '')

def findLibraryPath():
    """ Find the library path for HDU headers """
    path = os.path.split(os.path.abspath(__file__))[0]

    if os.path.exists(os.path.join(path, 'lib/header_primaryHDU.txt')):
        return os.path.join(path, 'lib')
    elif os.path.exists(os.path.join(path, 'header_primaryHDU.txt')):
        return path
    elif os.path.exists('header_primaryHDU.txt'):
        return './'
    else:
        raise IOError("Cannot find header files. Called from findLibraryPath() in sdfits.py")


def extractMid(x):
    """ Extract the mid part of an array """
    return x[len(x)/4:3*len(x)/4]

def fitLine(x, y, n_chans):
    """ Fit a line to data, only using central channels """
    x_fine = np.linspace(x[0], x[-1], n_chans)
    x = x[len(x)/4:3*len(x)/4]
    y = y[len(y)/4:3*len(y)/4]
    p = np.polyfit(x[::-1], y, 1)  # linear fit
    v = np.polyval(p, x_fine)
    return v
    
def loadDiodeTemp(h6, filename):
    """ Load a diode temp csv """
    
    f_fine = h6.freqs
    f      = h6.freqs_cal
    
    #temps_x = np.fromfile(filename_x).reshape([13,16])
    #temps_y = np.fromfile(filename_y).reshape([13,16])

    if filename.endswith('.hdf') or filename.endswith('h5') or filename.endswith('.hdf5'):
        temps = mbcal(filename)
    else:
        temps = np.fromfile(filename).reshape([26,16])
    temps_x = temps[0:13]
    temps_y = temps[13:26]

    temps_fine_x = np.zeros([13, 8192])
    temps_fine_y = np.zeros([13, 8192])
    
    for i in range(0,13):
        temps_fine_x[i] = fitLine(f, temps_x[i], 8192)
        temps_fine_y[i] = fitLine(f, temps_y[i], 8192)
        
    return temps_fine_x, temps_fine_y

def applyCal(beam, row, freqs, freqs_cal, cf, T_d_x, T_d_y):
    """ Apply basic calibration in Jy
    
    P_sys / (CF*P_on- CF* p_off) * T_d  = T_sys 
    
    P_sys: Total power in channel
    P_on:  Noise diode on
    P_off: Noise diode off
    T_d:   Temperature of the diode
    CF:    Cal factor which relates diode measurement P_on to P_sys
    
    """
    
    P_sys_xx = beam.cols.xx[row].astype('float')
    xx_on    = beam.cols.xx_cal_on[row].astype('float')
    xx_off   = beam.cols.xx_cal_off[row].astype('float')
    P_on_xx  = np.average(extractMid(xx_on))
    P_off_xx = np.average(extractMid(xx_off))
    
    #P_on_xx  = fitLine(freqs_cal, xx_on, len(freqs))
    #P_off_xx = fitLine(freqs_cal, xx_off, len(freqs))

    P_sys_yy = beam.cols.yy[row].astype('float')
    yy_on    = beam.cols.yy_cal_on[row].astype('float')
    yy_off   = beam.cols.yy_cal_off[row].astype('float')
    P_on_yy  = np.average(extractMid(yy_on))
    P_off_yy = np.average(extractMid(yy_off))
    
    #P_on_yy  = fitLine(freqs_cal, yy_on, len(freqs))
    #P_off_yy = fitLine(freqs_cal, yy_off, len(freqs))

    
    T_sys_xx = P_sys_xx / (cf*P_on_xx - cf*P_off_xx) * T_d_x
    T_sys_yy = P_sys_yy / (cf*P_on_yy - cf*P_off_yy) * T_d_y
    
    return T_sys_xx, T_sys_yy

def computeTsys(beam, row, T_d_x, T_d_y):
    """ Compute Tsys from frequency data 
    
    T_sys = Td / (Pon/Poff - 1)
    
    """
    
    xx_on    = beam.cols.xx_cal_on[row].astype('float')
    xx_off   = beam.cols.xx_cal_off[row].astype('float')
    
    yy_on    = beam.cols.yy_cal_on[row].astype('float')
    yy_off   = beam.cols.yy_cal_off[row].astype('float')
    
    T_sys_x = np.average(T_d_x[len(T_d_x)/4:3*len(T_d_x)/4]) / (xx_on/xx_off -1)
    T_sys_y = np.average(T_d_y[len(T_d_x)/4:3*len(T_d_x)/4]) / (yy_on/yy_off -1)
    
    l = len(T_sys_x)
    return np.average(T_sys_x[l/4:3*l/4]), np.average(T_sys_y[l/4:3*l/4])

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


def generateBlankSDFits(num_rows, header_primary, header_tbl, coldef_file):
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
    print header_primary
    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateBlankDataHDU(num_rows, header_tbl, coldef_file)
    
    # Insert creation date
    time_tuple = time.gmtime()
    date_str = time.strftime("%Y-%m-%dT%H:%M:%S", time_tuple)
    prhdu.header['DATE'] = date_str
        
    hdulist = pf.HDUList([prhdu, tbhdu])
    
    return hdulist

def generateSDFitsFromMbcorr(input_file, header_primary=None, header_tbl=None, coldef_file=None):
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

    if header_primary is None:
        header_primary = path + '/header_primaryHDU.txt'
    if header_tbl is None:
        header_tbl = path + '/header_dataHDU.txt'
    if coldef_file is None:
        coldef_file = path + '/coldefs_dataHDU.txt'

    prhdu = generatePrimaryHDU(header_primary)
    tbhdu = generateDataHDU(input_file, header_tbl, coldef_file)
    hdulist = pf.HDUList([prhdu, tbhdu])

    return hdulist

def generateSDFitsFromHipsr(filename_in, path_in, filename_out, path_out, write_stokes=0, cal=None):
    """ Generate an SD-FITS file from a hipsr5 file

    Parameters
    ----------
    filename_in: str
        Name of HDF5 file to load (input data)
    path_in: str
        path to HDF5 file (input data)
    filename_out: str
        Name for output SD-FITS file (output data)
    path_out: str
        Path to output SD-FITS file (output data)
    write_stokes: 0, 1, or 2 (int)
        What kind of data to write to file
        0 write autocorrs only,
        1 write cross-correlations too, stored in XPOLDATA
        2 write Stokes I Q U V
    cal: None or str
        Path to calibration file. If None, will use defaults.
    """
    
    # Open h5 file
    print "\nOpening files"
    print "-------------"
    h5file = os.path.join(path_in, filename_in)
    out_file = os.path.join(path_out, filename_out)
    h6 = Hipsr6(h5file)
    pointing = h6.tb_pointing.cols
    obs      = h6.tb_observation.cols
    obs_mode = obs.obs_mode[0].strip()
    ref_beams= obs.ref_beam[:]

    freqs     = h6.freqs
    freqs_cal = h6.freqs_cal
    
    print "Input file: %s"%h6.h5.filename
    print h6

    if cal == None:
        abspath = os.path.abspath( __file__ ).replace('sdfits.pyc', '').replace('sdfits.py', '')
        #diode_cal_file_x  = "%s/diode_jy_x.cal"%abspath
        #diode_cal_file_y  = "%s/diode_jy_y.cal"%abspath
        diode_cal_file = "%s/diode_jy.cal"%abspath
    else:
        diode_cal_file = cal

    print "Using calibration %s"%cal
    diode_temps_x, diode_temps_y = loadDiodeTemp(h6, diode_cal_file)

    scan_pointing_len = h6.tb_scan_pointing.shape[0]
    
    tb_lengths = []
    for beam in h6.h5.root.raw_data:
        if beam.shape[0] != scan_pointing_len:
            beam_id = int(beam.name.lstrip('beam_'))
            print "WARNING: beam %i len: %i, scan_pointing len: %i"%(beam_id, beam.shape[0], scan_pointing_len)
        tb_lengths.append(np.min([beam.shape[0], scan_pointing_len]))
        
     
    num_acc  = np.max(tb_lengths) 
    num_rows = num_acc * 13



    if num_acc == 0:
        print "No data in %s. Skipping."%h5file
        return -1
    
    print "No accumulations: %s, no rows: %s"%(num_acc, num_rows)

    # We now need to generate a blank SD-FITS file, with the same number of rows
    print "\nGenerating blank SD-FITS file with %i rows..."%num_rows

    path = findLibraryPath()
    if obs_mode == 'MXCAL':
        header_primary = os.path.join(path, 'header_primaryHDU.txt')
        header_tbl = os.path.join(path, 'header_dataHDU_mxcal.txt')
        coldef_file = os.path.join(path, 'coldefs_dataHDU_mxcal.txt')
    elif write_stokes == 2:
        print "Stokes flag found - writing I,Q,U,V"
        header_primary = os.path.join(path, 'header_primaryHDU.txt')
        header_tbl = os.path.join(path, 'header_dataHDU_stokes.txt')
        coldef_file = os.path.join(path, 'coldefs_dataHDU_stokes.txt')
    elif write_stokes == 0:
        print "Writing XX, YY"
        header_primary = os.path.join(path, 'header_primaryHDU.txt')
        header_tbl = os.path.join(path, 'header_dataHDU.txt')
        coldef_file = os.path.join(path, 'coldefs_dataHDU.txt')
    else:
        print "Writing XX, YY, XY, YX"
        header_primary = os.path.join(path, 'header_primaryHDU.txt')
        header_tbl = os.path.join(path, 'header_dataHDU_xpol.txt')
        coldef_file = os.path.join(path, 'coldefs_dataHDU_xpol.txt')
    
    hdulist = generateBlankSDFits(num_rows, header_primary, header_tbl, coldef_file)
    print hdulist.info()
    
    # Next, we copy over observation data    
    print "Filling new SD-FITS with HIPSR data..."
    sdtab    = hdulist[1].data
    sdhead   = hdulist[1].header

    # Fill in header values
    sdhead["OBSERVER"] = obs.observer[0]
    sdhead["PROJID"]   = obs.project_id[0]
    
    # Fill in common values
    # NEW METHOD OF TIMESTAMPING - AUG 27 2013
    ref_time  = int(h6.h5.root.raw_data.beam_01.cols.timestamp[0])
    ref_id    = int(h6.h5.root.raw_data.beam_01.cols.id[0])
    ref_clk   = 800e6 # Clock frequency 800 MHz
    num_chans = h6.h5.root.raw_data.beam_01.cols.xx[0].shape[0]
    acc_len   = h6.h5.root.firmware_config.cols.acc_len[0]
    ref_delta = num_chans * acc_len * 2 / ref_clk

    print "Filling in common values... ",
    sdtab["SCAN"][:]     = 1
    sdtab["EXPOSURE"][:] = ref_delta
    sdtab["OBJECT"][:]   = pointing.source[0]
    sdtab["OBJ-RA"][:]   = pointing.ra[0]
    sdtab["OBJ-DEC"][:]  = pointing.dec[0]
    sdtab["RESTFRQ"][:]  = obs.frequency[0] * 1e6
    sdtab["FREQRES"][:]  = np.abs(obs.bandwidth[0])*1e6 / num_chans
    sdtab["BANDWID"][:]  = np.abs(obs.bandwidth[0]) * 1e6
    sdtab["CRPIX1"][:]   = num_chans/2 + 1
    sdtab["CRVAL1"][:]   = obs.frequency[0] * 1e6
    sdtab["CDELT1"][:]   = np.abs(obs.bandwidth[0])*1e6 / num_chans
    sdtab["FLAGGED"][:]  = 0
    sdtab["SCANRATE"][:] = obs.scan_rate[0] / 60 # Deg/min to deg/s


    # TCS INFO
    sdtab["OBSMODE"][:]  = obs.obs_mode[0] 
    sdtab["IF"][:]       = 1
    print "OK."
    
    row_sd   = 0
    cycle_id = 0
    
    flipped = False
    if obs.bandwidth[0] < 0:
        flipped = True
    
    print "Filling in unique values... "
    num_cycles = np.min([scan_pointing_len, num_acc])
    for row_h5 in range(num_acc):
        cycle_id += 1 # Starts at 1 in SD-FITS file

        for beam in h6.h5.root.raw_data:
            beam_id = int(beam.name.lstrip('beam_'))
            LinePrint("%i of %i"%(row_sd, num_rows))
            
            if cycle_id <= num_cycles:
                raj_id = "mb%s_raj"%beam.name.lstrip('beam_')
                dcj_id = "mb%s_dcj"%beam.name.lstrip('beam_')
                
                sdtab["CYCLE"][row_sd]   = cycle_id

                # Fix beam mapping (remove after fixing mapping)
                sdtab["BEAM"][row_sd]     = beam_id
                
                sdtab["CRVAL3"][row_sd]   = h6.tb_scan_pointing.col(raj_id)[cycle_id-1]
                sdtab["CRVAL4"][row_sd]   = h6.tb_scan_pointing.col(dcj_id)[cycle_id-1]

                # AZ, EL and PARANGLE should be stored for beam 1 only
                if beam_id == 1:
                    sdtab["AZIMUTH"][row_sd]  = h6.tb_scan_pointing.col("azimuth")[cycle_id-1]
                    sdtab["ELEVATIO"][row_sd] = h6.tb_scan_pointing.col("elevation")[cycle_id-1]
                    sdtab["PARANGLE"][row_sd] = h6.tb_scan_pointing.col("par_angle")[cycle_id-1]

                #sdtab["FOCUSAXI"][row_sd] = h6.tb_scan_pointing.col("focus_axi")[cycle_id-1]
                sdtab["FOCUSTAN"][row_sd] = h6.tb_scan_pointing.col("focus_tan")[cycle_id-1]

                # This is confusing - but it looks like FOCUSROT should be 15.0, which is sent as feed_angle
                # Likewise, focusaxi is probably supposed to be what we receive as focus_rot
                focus_rot = h6.tb_scan_pointing.col("focus_rot")[cycle_id-1]
                sdtab["FOCUSROT"][row_sd] = focus_rot
                sdtab["FOCUSAXI"][row_sd] = h6.tb_observation.col("feed_angle")[0]

                try:

                    # OLD - 27 Aug 2013
                    #timestamp  = beam.cols.timestamp[row_h5]
                    # New - based off integration length
                    if beam_id == 1:
                        new_id = beam.cols.id[row_h5]
                        timestamp = (new_id - ref_id) * ref_delta + ref_time
                        date_obs, time = timestamp2dt(timestamp)

                    sdtab["DATE-OBS"][row_sd] = date_obs
                    sdtab["TIME"][row_sd]     = time

                    ref_beam = np.min(np.abs(timestamp - ref_beams[beam_id-1]))
                    
                    # Compute T_sys for each beam
                    T_d_x = diode_temps_x[beam_id-1]
                    T_d_y = diode_temps_y[beam_id-1]
                    T_sys_x, T_sys_y = computeTsys(beam, row_h5, T_d_x, T_d_y)
                
                    #print T_sys_x, T_sys_y
                    sdtab["TSYS"][row_sd] = (T_sys_x, T_sys_y)
                    sdtab["TCAL"][row_sd] = (np.average(extractMid(T_d_x)), np.average(extractMid(T_d_y)))
                    #sdtab["CALFCTR"][row_sd] = (1, 1)

                    if obs_mode == 'MXCAL':
                        sdtab["REFBEAM"][row_sd] = ref_beam

                    if write_stokes == 2:
                        # Currently not calibrating!
                        xx = beam.cols.xx[row_h5].astype('float32') 
                        yy = beam.cols.yy[row_h5].astype('float32') 
                        re_xy = beam.cols.re_xy[row_h5].astype('float32') 
                        im_xy = beam.cols.im_xy[row_h5].astype('float32')
                        
                        # Blank DC bin
                        xx[0], yy[0],re_xy[0], im_xy[0] = np.zeros(4)
                        if flipped:
                            xx, yy, re_xy, im_xy = xx[::-1], yy[::-1], re_xy[::-1], im_xy[::-1]
                        
                        
                        xx = xx / np.average(extractMid(xx)) * T_sys_x
                        yy = yy / np.average(extractMid(yy)) * T_sys_y
                        re_xy = re_xy / np.average(extractMid(re_xy)) * np.sqrt(T_sys_x * T_sys_y)
                        im_xy = im_xy / np.average(extractMid(im_xy)) * np.sqrt(T_sys_x * T_sys_y)
                        
                        # Ettore tells me Parkes uses this definition
                        # i.e. that I is the average of xx + yy
                        ii = (xx + yy) / 2
                        qq = (xx - yy) / 2
                        uu = re_xy
                        vv = im_xy
                        
                        # Form one data vector
                        data1 = np.append(ii, qq)
                        data2 = np.append(uu, vv)
                        data  = np.append(data1, data2)
                        data  = data.reshape([1,1,4,num_chans])
                    else:
                        
                        xx = beam.cols.xx[row_h5].astype('float32') 
                        yy = beam.cols.yy[row_h5].astype('float32') 
                        # Blank DC bin
                        xx[0], yy[0] = 0,0

                        if write_stokes == 1:
                            re_xy = beam.cols.re_xy[row_h5].astype('float32')
                            im_xy = beam.cols.im_xy[row_h5].astype('float32')
                            re_xy = re_xy / np.average(extractMid(re_xy)) * np.sqrt(T_sys_x * T_sys_y)
                            im_xy = im_xy / np.average(extractMid(im_xy)) * np.sqrt(T_sys_x * T_sys_y)
                            re_xy[0], im_xy[0] = 0, 0

                        if flipped:
                            xx, yy = xx[::-1], yy[::-1]                           
                            if write_stokes:
                                re_xy, im_xy = re_xy[::-1], im_xy[::-1]

                        #print "cal factor: %2.3f"%cf
                        #print "Diode temp: %s"%T_d
                        #xx, yy = applyCal(beam, row_h5, freqs, freqs_cal, cf, T_d_x, T_d_y)
                        
                        xx = xx / np.average(extractMid(xx)) * T_sys_x
                        yy = yy / np.average(extractMid(yy)) * T_sys_y
                        
                        # Multibeam stats screws up if it encounters division by 1
                        xx[xx <= 1 ] = 1  
                        yy[yy <= 1 ] = 1 
                        
                        do_flagger = True
                        if do_flagger:
                            flags = np.zeros(len(xx))
                            #flags[xx>T_sys_x*5] = 1
                            #flags[yy>T_sys_x*5] = 1
                            flags[xx==1] = 1
                            flags[yy==1] = 1
                            flags = np.append(flags, flags)
                            flags = flags.reshape([1,1,2,num_chans])
                            sdtab["FLAGGED"][row_sd] = flags
                        
                        data = np.append(xx, yy)
                        data = data.reshape([1,1,2,num_chans])
                    
                    sdtab["DATA"][row_sd] = data

                    if write_stokes == 1:
                        sdtab["XPOLDATA"][row_sd] = np.row_stack((re_xy, im_xy)).flatten()
                    
                except:
                    if beam.name != 'beam_02':
                        print "\nWARNING: missing row in %s"%beam.name
                        print "Current index: %i"%row_h5
                        print "Row length: %i"%beam.shape[0]
                        raise
                    try:
                        sdtab["FLAGGED"][row_sd] = np.ones_like([1,1,2,num_chans])
                    except ValueError:
                        pass
                row_sd += 1
            else:
                print "WARNING: scan_pointing table is not complete."
                print "%s table length: %i"%(beam.name, beam.shape[0])
                print "scan_pointing table length: %i"%scan_pointing_len

    
    h6.h5.close()
    
    if os.path.exists(out_file):
        print "\nInfo: File exists, deleting..."
        os.remove(out_file)

    print "\nInfo: Saving to file"
    hdulist.writeto(out_file)
    hdulist.close()

