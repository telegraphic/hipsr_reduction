#!/usr/bin/env python
"""
hipsr-converter.py
==================

This script starts a graphical user interface for converting HIPSR + MBCORR data to SD-FITS.
Use this script when converting HIPSR + MBCORR data taken in tandem.

"""

# Imports
import sys
from lib.sdfits import *
from lib import eager_weaver

try:
    from termcolor import cprint
except ImportError:
    def cprint_fallback(textstr, color):
        print textstr

    cprint = cprint_fallback


# Python metadata
__version__  = "v2.0 - Ballistic Bandicoot"
__author__   = "Danny Price"
__email__    = "dprice@cfa.harvard.edu"
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))

try:
    import lib.qt_compat as qt_compat
    QtGui = qt_compat.import_module("QtGui")
    QtCore = qt_compat.QtCore
    USES_PYSIDE = qt_compat.is_pyside()

except:
    print "Error: cannot load PySide or PyQt4. Please check your install."
    exit()

try:
    import numpy as np
except:
    print "Error: cannot load Numpy. Please check your install."
    exit()

try:
    import pyfits as pf
except ImportError:
    try:
        from astropy.io import fits as pf
        print "Using Astropy for FITS I/O"
    except:
        print "Error: cannot load PyFITS or AstroPY I/O. Please check your install."
        exit()

try:
    import tables as tb
except:
    print "Error: cannot load PyTables. Please check your install."
    exit()


class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        last_in, last_mb, last_out = self.load_last()

        self.in_combox = self.createComboBox(last_in)
        self.in_label  = QtGui.QLabel("HIPSR input directory:")
        self.in_browse = self.createButton("&Browse...", self.in_set)
        self.in_label.setToolTip("Select input directory (HDF files)")
        self.in_combox.setToolTip("Select input directory (HDF files)")

        self.mb_combox = self.createComboBox(last_mb)
        self.mb_label  = QtGui.QLabel("MBCORR input directory:")
        self.mb_browse = self.createButton("&Browse...", self.mb_set)
        self.mb_label.setToolTip("Select MBCORR input directory (SD-FITS files)")
        self.mb_combox.setToolTip("Select MBCORR input directory (SD-FITS files)")


        self.out_combox = self.createComboBox(last_out)
        self.out_label  = QtGui.QLabel("Output directory:")
        self.out_browse = self.createButton("&Browse...", self.out_set)
        self.out_label.setToolTip("Select output directory (SD-FITS)")
        self.out_combox.setToolTip("Select output directory (SD-FITS")
        
        self.convert_button = self.createButton("&Convert", self.convert)
        
        #self.rb_autos = QtGui.QRadioButton("Write autocorrs", self)
        #self.rb_xpol    = QtGui.QRadioButton("Write cross-pol", self)
        #self.rb_stokes  = QtGui.QRadioButton("Write Stokes", self)

        #self.rb_autos.setChecked(True)

        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.in_label, 0, 0)
        mainLayout.addWidget(self.in_combox, 0, 1)
        mainLayout.addWidget(self.in_browse, 0, 2)
        mainLayout.addWidget(self.mb_label, 1, 0)
        mainLayout.addWidget(self.mb_combox, 1, 1)
        mainLayout.addWidget(self.mb_browse, 1, 2)
        mainLayout.addWidget(self.out_label, 2, 0)
        mainLayout.addWidget(self.out_combox, 2, 1)
        mainLayout.addWidget(self.out_browse, 2, 2)
        #mainLayout.addWidget(self.rb_autos, 3, 1)
        #mainLayout.addWidget(self.rb_xpol, 4, 1)
        #mainLayout.addWidget(self.rb_stokes, 5, 1)
        mainLayout.addWidget(self.convert_button, 3, 2)

        self.setLayout(mainLayout)

        self.setWindowTitle("HIPSR-MBCORR tandem observation data converter")

    def load_last(self):
        try:
            f = open(QtCore.QDir.currentPath()+'/.last_tandem')
            last_in = f.readline().strip('\n')
            last_mb = f.readline().strip('\n')
            last_out = f.readline().strip('\n')
            f.close()
            if os.path.exists(last_in) and os.path.exists(last_out):
                return last_in, last_mb, last_out
            else:
                raise IOError
        except:
            return QtCore.QDir.currentPath(), QtCore.QDir.currentPath(), QtCore.QDir.currentPath()

    def save_last(self):
        try:
            f = open(QtCore.QDir.currentPath()+'/.last_tandem', 'w')
            f.write(self.in_combox.currentText()+'\n')
            f.write(self.mb_combox.currentText()+'\n')
            f.write(self.out_combox.currentText()+'\n')
            f.close()
        except IOError:
            pass

    def in_set(self):
        last_in, last_mb, last_out = self.load_last()
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select HIPSR input directory",
                                                           last_in + '/..')
        if directory:
            if self.in_combox.findText(directory) == -1:
                self.in_combox.addItem(directory)
            self.in_combox.setCurrentIndex(self.in_combox.findText(directory))

    def mb_set(self):
        last_in, last_mb, last_out = self.load_last()
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select MBCORR input directory",
                                                           last_mb + '/..')
        if directory:
            if self.mb_combox.findText(directory) == -1:
                self.mb_combox.addItem(directory)
            self.mb_combox.setCurrentIndex(self.mb_combox.findText(directory))


    def out_set(self):
        last_in, last_mb, last_out = self.load_last()
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select SD-FITS ouput directory",
                                                           last_out + '/..')
        if directory:
            if self.out_combox.findText(directory) == -1:
                self.out_combox.addItem(directory)
            self.out_combox.setCurrentIndex(self.out_combox.findText(directory))    

    def updateComboBox(comboBox):
        if comboBox.findText(comboBox.currentText()) == -1:
            comboBox.addItem(comboBox.currentText())

    def createButton(self, text, member):
        button = QtGui.QPushButton(text)
        button.clicked.connect(member)
        return button

    def createComboBox(self, text=""):
        comboBox = QtGui.QComboBox()
        comboBox.setEditable(True)
        comboBox.addItem(text)
        comboBox.setSizePolicy(QtGui.QSizePolicy.Expanding,
                QtGui.QSizePolicy.Preferred)
        return comboBox

    def convert(self):

        self.save_last()

        print("HIPSR-MBCORR tandem converter")
        print("-----------------------------")
        print("Input directory (HIPSR): %s"%self.in_combox.currentText())
        print("Input directory (MBCORR): %s"%self.mb_combox.currentText())
        print("Output directory: %s"%self.out_combox.currentText())

        hipsr_dir  = self.in_combox.currentText()
        mbcorr_dir = self.mb_combox.currentText()
        mbcorr_files = eager_weaver.findMbcorrFiles(self.mb_combox.currentText())
        output_dir   = self.out_combox.currentText()

        # Make sure output directory exists
        if not os.path.exists(output_dir):
            print("Creating directory %s"%output_dir)
            os.makedirs(output_dir)

        i = 0
        for mb_filename in mbcorr_files:
            i += 1
            cprint("\nfile %i of %i (%02d%%)"%(i, len(mbcorr_files), float(i)/len(mbcorr_files)*100), 'green')
            cprint("-------------------", 'green')

            hp_filename, t_diff = eager_weaver.filePairer(mb_filename, hipsr_dir)
            if t_diff >= 60:
               print "No match found for %s"%mb_filename
               break

            if t_diff <= 60:
                print "MBCORR input file:     %s"%mb_filename
                print "Closest matching file: %s"%hp_filename
                print "Time delta: %d\n"%t_diff

                out_filename = os.path.join(output_dir, 'hipsr_'+os.path.basename(mb_filename))
                eager_weaver.eagerWeaver(mb_filename, hp_filename, out_filename,
                                         hp_search_dir=hipsr_dir, sd_search_dir=mbcorr_dir, gmt_diff=0)
            else:
                print "No matching file found. Skipping..."



        print("DONE!")


if __name__ == '__main__':

    import sys

    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()