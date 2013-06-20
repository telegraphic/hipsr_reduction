#!/usr/bin/env python
"""
hipsr-converter.py
==================

This script starts a graphical user interface for converting HIPSR data to SD-FITS.

"""



# Imports
import sys
from optparse import OptionParser

import hipsr_core.config as config
from hipsr_core.sdfits import *

# Python metadata
__version__  = config.__version__
__author__   = config.__author__
__email__    = config.__email__
__license__  = config.__license__
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))


try:
    import hipsr_core.qt_compat as qt_compat
    QtGui = qt_compat.import_module("QtGui")
    QtCore = qt_compat.QtCore
    
    USES_PYSIDE = qt_compat.is_pyside()
    
    #import PyQt4
    #from PyQt4 import QtGui, QtCore
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
except:
    print "Error: cannot load PyFITS. Please check your install."
    exit()

try:    
    import tables as tb
except:
    print "Error: cannot load PyTables. Please check your install."
    exit()

class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        
        self.in_combox = self.createComboBox(QtCore.QDir.currentPath())
        self.in_label  = QtGui.QLabel("Input directory:")
        self.in_browse = self.createButton("&Browse...", self.in_set)
        self.in_label.setToolTip("Select input directory (HDF files)")
        self.in_combox.setToolTip("Select input directory (HDF files)")

        self.out_combox = self.createComboBox(QtCore.QDir.currentPath())
        self.out_label  = QtGui.QLabel("Output directory:")
        self.out_browse = self.createButton("&Browse...", self.out_set)
        self.out_label.setToolTip("Select output directory (SD-FITS)")
        self.out_combox.setToolTip("Select output directory (SD-FITS")
        
        self.convert_button = self.createButton("&Convert", self.convert)
        
        self.cb_stokes = QtGui.QCheckBox('Write Stokes?', self)
        
        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.in_label, 0, 0)
        mainLayout.addWidget(self.in_combox, 0, 1)
        mainLayout.addWidget(self.in_browse, 0, 2)
        mainLayout.addWidget(self.out_label, 1, 0)
        mainLayout.addWidget(self.out_combox, 1, 1)
        mainLayout.addWidget(self.out_browse, 1, 2)
        mainLayout.addWidget(self.cb_stokes, 2, 1)
        mainLayout.addWidget(self.convert_button, 2, 2)
        
        self.setLayout(mainLayout)

        self.setWindowTitle("HIPSR-converter: HDF5 to SD-FITS")
    

    def in_set(self):
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select HDF input directory",
                QtCore.QDir.currentPath())

        if directory:
            if self.in_combox.findText(directory) == -1:
                self.in_combox.addItem(directory)

            self.in_combox.setCurrentIndex(self.in_combox.findText(directory))

    def out_set(self):
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select HDF input directory",
                QtCore.QDir.currentPath())

        if directory:
            if self.out_combox.findText(directory) == -1:
                self.out_combox.addItem(directory)

            self.out_combox.setCurrentIndex(self.out_combox.findText(directory))    

    @staticmethod
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
        
        print "HIPSR SD-FITS writer"
        print "--------------------"
        print "Input directory: %s"%self.in_combox.currentText()
        print "Output directory: %s"%self.out_combox.currentText()
    
        # Regular expression to match h5 extension
        regex = '([0-9A-Za-z-_]+).hdf'
        filelist = []
        path = self.in_combox.currentText()
        out_path = self.out_combox.currentText()
        for filename in os.listdir(path):
            match = re.search(regex, filename)
            if match:
                filelist.append(filename)
    
        # Make sure output directory exists
        if not os.path.exists(out_path):
            print "Creating directory %s"%out_path
            os.makedirs(out_path)
    
        i = 1
        
        if self.cb_stokes.isChecked(): ws = True
        else: ws = False
            
        for file_in in filelist:
            print "Creating file %i of %i... \n"%(i, len(filelist))
            file_out = file_in.rstrip('.h5').rstrip('.hdf') + '.sdfits'
                
            generateSDFitsFromHipsr(file_in, path, file_out, out_path, write_stokes=ws)
        
            i += 1
    
        print "DONE!"


if __name__ == '__main__':

    import sys

    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
