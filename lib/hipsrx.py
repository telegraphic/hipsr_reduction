import numpy as np
import tables as tb
from datetime import datetime

class Hipsr6(object):
    """ Hipsr6 data object thingo """
    def __init__(self, h5_filename):
        
        self.h5 = tb.openFile(h5_filename)
        
        # Data access helpers
        self.beams              = [b for b in self.h5.root.raw_data]
        self.tb_observation     = self.h5.root.observation
        self.tb_pointing        = self.h5.root.pointing
        self.tb_scan_pointing   = self.h5.root.scan_pointing
        self.tb_firmware_config = self.h5.root.firmware_config
        self.tb_weather         = self.h5.root.weather
        
        # Things that need to be computed
        self.computeFreqs()
    
    def close(self):
        """ Close hdf file """
        self.h5.close()
        
    def computeFreqs(self):
        """ Compute frequency channels from file info"""
        # Compute channel frequencies
        bw          = self.h5.root.observation[0]["bandwidth"]
        cent_fr     = self.h5.root.observation[0]["frequency"]
        n_chans     = self.h5.root.raw_data.beam_01.coldtypes['xx'].shape[0]
        n_chans_cal = self.h5.root.raw_data.beam_01.coldtypes['xx_cal_on'].shape[0]     
        
        self.freqs     = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, n_chans)
        self.freqs_cal = np.linspace(cent_fr - np.abs(bw)/2, cent_fr + np.abs(bw)/2, n_chans_cal)
        
        self.flipped = True if bw < 0 else False 
        if self.flipped: self.freqs, self.freqs_cal = self.freqs[::-1], self.freqs_cal[::-1]
    
    def __repr__(self):
        to_print = 'HIPSR6 file info:\n'
        to_print += 'Filename:     %s\n'%self.h5.filename
        to_print += 'Timestamp:    %s\n'%datetime.fromtimestamp(self.h5.root.pointing.cols.timestamp[0])
        to_print += 'First object: %s\n'%self.h5.root.pointing.cols.source[0]
        return to_print