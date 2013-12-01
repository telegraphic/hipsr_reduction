#!/usr/bin/env python
"""
hipsr-converter.py
==================

This script starts a graphical user interface for converting HIPSR data to SD-FITS.

"""

# Imports
import sys
from lib.sdfits import *

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

        last_in, last_out = self.load_last()

        self.in_combox = self.createComboBox(last_in)
        self.in_label  = QtGui.QLabel("Input directory:")
        self.in_browse = self.createButton("&Browse...", self.in_set)
        self.in_label.setToolTip("Select input directory (HDF files)")
        self.in_combox.setToolTip("Select input directory (HDF files)")

        self.out_combox = self.createComboBox(last_out)
        self.out_label  = QtGui.QLabel("Output directory:")
        self.out_browse = self.createButton("&Browse...", self.out_set)
        self.out_label.setToolTip("Select output directory (SD-FITS)")
        self.out_combox.setToolTip("Select output directory (SD-FITS")
        
        self.convert_button = self.createButton("&Convert", self.convert)
        
        self.rb_autos = QtGui.QRadioButton("Write autocorrs", self)
        self.rb_xpol    = QtGui.QRadioButton("Write cross-pol", self)
        self.rb_stokes  = QtGui.QRadioButton("Write Stokes", self)

        self.rb_autos.setChecked(True)

        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.in_label, 0, 0)
        mainLayout.addWidget(self.in_combox, 0, 1)
        mainLayout.addWidget(self.in_browse, 0, 2)
        mainLayout.addWidget(self.out_label, 1, 0)
        mainLayout.addWidget(self.out_combox, 1, 1)
        mainLayout.addWidget(self.out_browse, 1, 2)
        mainLayout.addWidget(self.rb_autos, 2, 1)
        mainLayout.addWidget(self.rb_xpol, 3, 1)
        mainLayout.addWidget(self.rb_stokes, 4, 1)
        mainLayout.addWidget(self.convert_button, 5, 2)

        self.setLayout(mainLayout)

        self.setWindowTitle("HIPSR-converter: HDF5 to SD-FITS")

    def load_last(self):
        try:
            f = open(QtCore.QDir.currentPath()+'/.last')
            last_in = f.readline().strip('\n')
            last_out = f.readline().strip('\n')
            f.close()
            if os.path.exists(last_in) and os.path.exists(last_out):
                return last_in, last_out
            else:
                raise IOError
        except:
            return QtCore.QDir.currentPath(), QtCore.QDir.currentPath()

    def save_last(self):
        try:
            f = open(QtCore.QDir.currentPath()+'/.last', 'w')
            f.write(self.in_combox.currentText()+'\n')
            f.write(self.out_combox.currentText()+'\n')
            f.close()
        except IOError:
            pass

    def in_set(self):
        last_in, last_out = self.load_last()
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select HDF input directory",
                                                           last_in + '/..')

        if directory:
            if self.in_combox.findText(directory) == -1:
                self.in_combox.addItem(directory)

            self.in_combox.setCurrentIndex(self.in_combox.findText(directory))

    def out_set(self):
        last_in, last_out = self.load_last()
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select HDF input directory",
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

        print("HIPSR SD-FITS writer")
        print("--------------------")
        print("Input directory: %s"%self.in_combox.currentText())
        print("Output directory: %s"%self.out_combox.currentText())

        # Regular expression to match h5 extension
        regex = '([0-9A-Za-z-_]+).(hdf|h5)'
        filelist = []
        path = self.in_combox.currentText()
        out_path = self.out_combox.currentText()
        for filename in os.listdir(path):
            match = re.search(regex, filename)
            if match:
                filelist.append(filename)

        # Make sure output directory exists
        if not os.path.exists(out_path):
            print("Creating directory %s"%out_path)
            os.makedirs(out_path)

        i = 1

        # Check what type of data is to be written
        ws = 0
        if self.rb_xpol.isChecked():
            ws = 1
        if self.rb_stokes.isChecked():
            ws = 2

        for file_in in filelist:
            print("Creating file %i of %i... \n"%(i, len(filelist)))
            file_out = os.path.splitext(file_in)[0] + '.sdfits'
            time.sleep(1)
            generateSDFitsFromHipsr(file_in, path, file_out, out_path, write_stokes=ws)

            i += 1

        print("DONE!")


if __name__ == '__main__':

    import sys

    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()