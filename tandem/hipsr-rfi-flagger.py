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
    import lib.qt_compat as qt_compat
    QtCore = qt_compat.QtCore
    QtGui = qt_compat.import_module("QtGui")
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

from lib.printers import LinePrint

def rfiFlagger(sd_filename):
    """ Reverse and flag a file """
    #print "opening %s"%sd_filename
    sd = pf.open(sd_filename, mode='update')
    sd_len = sd[1].data.shape[0]

    i = 0
    for row in range(sd_len):
        i += 1
        x = sd[1].data["DATA"][0,0,0,0]
        y = sd[1].data["DATA"][0,0,0,1]
        f = sd[1].data["FLAGGED"][0,0,0,0]
        d = np.abs(x - y)
        avg = np.average(d[f==0])
        std = np.std(d[f==0])

        f2 = d < avg + 3*std

        fnew = f | f2

        sd[1].data["FLAGGED"][row,0,0,0] = fnew.astype('int8')  #1345 MHz

    sd.flush()
    sd.close()

class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        
        self.in_combox = self.createComboBox('/home/dprice/data/P641/hipsr_flagged_sdfits')
        self.in_label  = QtGui.QLabel("SD-FITS directory:")
        self.in_browse = self.createButton("&Browse...", self.in_set)
        self.in_label.setToolTip("Select input directory (HDF files)")
        self.in_combox.setToolTip("Select input directory (HDF files)")

        self.convert_button = self.createButton("&Flag RFI", self.convert)

        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.in_label, 0, 0)
        mainLayout.addWidget(self.in_combox, 0, 1)
        mainLayout.addWidget(self.in_browse, 0, 2)
        mainLayout.addWidget(self.convert_button, 2, 1)
        self.setLayout(mainLayout)

        self.setWindowTitle("HIPSR-flipper")
    

    def in_set(self):
        directory = QtGui.QFileDialog.getExistingDirectory(self, "Select SD-FITS directory",
                QtCore.QDir.currentPath())

        if directory:
            if self.in_combox.findText(directory) == -1:
                self.in_combox.addItem(directory)

            self.in_combox.setCurrentIndex(self.in_combox.findText(directory))

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
        
        print "HIPSR SD-FITS flipper"
        print "--------------------"
        print "Input directory: %s"%self.in_combox.currentText()

        # Regular expression to match h5 extension
        regex = '([0-9A-Za-z-_]+).sdfits'
        filelist = []
        path = self.in_combox.currentText()
        for filename in os.listdir(path):
            match = re.search(regex, filename)
            if match:
                filelist.append(filename)

        i = 1
        for file_in in filelist:
            print "(%02d%%) Flagging %i of %i"%(float(i)/len(filelist)*100, i, len(filelist))
            rfiFlagger(os.path.join(self.in_combox.currentText(),file_in))
            i += 1
    
        print "DONE!"


if __name__ == '__main__':

    import sys

    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
