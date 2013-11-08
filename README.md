hipsr_reduction
===============

Data reduction scripts to convert HIPSR data into SD-FITS files that can be read by Livedata or similar
packages. The two main scripts are:

* hipsr-converter.py -- for conversion of HIPSR data into SD-FITS files.
* hipsr-converter-tandem.py -- for conversion of HIPSR + MBCORR data into combined SD-FITS files.


Requirements
------------

* numpy (>1.6) - older versions may cause issues!
* PyQt4 or PySide - Required for GUI. PyQt4 is easier to install on Linux, PySide is easier on Mac.
* pyfits (or astropy.io) - required for FITS reading/writing
* pytables - required for HDF5 reading / writing (HIPSR storage format)

Install notes
-------------

On Mac I've had success installing dependencies using homebrew (http://brew.sh) and pip (http://www.pip-installer.org).
On newish versions of Ubuntu, the apt-get sources should be new enough.

Usage
-----

To open the GUI, just type:

    python hipsr-converter.py
    
and select the input and output directories, then click 'convert'. The output of this script will
be sdfits files, that can then be run through Livedata / Gridzilla or similar to create science products.

Alternatively, if dealing with tandem MBCORR-HIPSR observations, run:

    python hipsr-converter-tandem.py

Note that the tandem converter requires the MBCORR directory to contain SD-FITS formatted files.
To get these, you'll have to run the *.rpf / *.mbf files through Livedata WITHOUT applying bandpass
calibration. For example, a step-by-step reduction would go something like:

1. Convert *.rpf files into *.sdfits with Livedata, without applying bandpass calibration.
2. Run hipsr-converter-tandem.py and stitch together sdfits files with HIPSR data (hdf)
3. Run the output files (which will be sdfits) through Livedata to apply bandpass calibration.
4. These calibrated sdfits files can then be loaded into gridzilla or similar to make a data cube.

Note that if you run a tandem observation and then try to calibrate without using MBCORR, you will get
unusable data. This is as MBCORR provides calibration information for tandem observation modes.
