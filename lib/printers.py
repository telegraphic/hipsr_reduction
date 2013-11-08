# encoding: utf-8
"""
lineprint.py
=============

Line printer - prints on same line instead of a new line each time.
""" 

__version__ = "0.2"

import sys, os, time
from datetime import datetime

class LinePrint():
    """
    Print things to stdout on one line dynamically
    """
    def __init__(self,data):
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()
        
class Logger(object):
    """ A logger which records the output of stdin to file """
    def __init__(self, filename, filepath):
        self.terminal = sys.stdout
        self.log      = open(os.path.join(filepath,filename), 'a')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message) 