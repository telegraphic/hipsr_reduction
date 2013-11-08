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
