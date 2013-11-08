#!/usr/bin/env python
"""
qt_compat.py
============

Compatibility script to fall back to PyQt4 if PySide isn't installed/working.

[1] http://qt-project.org/wiki/Differences_Between_PySide_and_PyQt
"""

from __future__ import with_statement
from __future__ import division

_TRY_PYSIDE = True

try:
	if not _TRY_PYSIDE:
		raise ImportError()
	import PySide.QtCore as _QtCore
	QtCore = _QtCore
	USES_PYSIDE = True
except ImportError:
	import sip
	sip.setapi('QString', 2)
	sip.setapi('QVariant', 2)
	import PyQt4.QtCore as _QtCore
	QtCore = _QtCore
	USES_PYSIDE = False


def _pyside_import_module(moduleName):
	pyside = __import__('PySide', globals(), locals(), [moduleName], -1)
	return getattr(pyside, moduleName)


def _pyqt4_import_module(moduleName):
	pyside = __import__('PyQt4', globals(), locals(), [moduleName], -1)
	return getattr(pyside, moduleName)

if USES_PYSIDE:
	import_module = _pyside_import_module

	Signal = QtCore.Signal
	Slot = QtCore.Slot
	Property = QtCore.Property
else:
	import_module = _pyqt4_import_module

	Signal = QtCore.pyqtSignal
	Slot = QtCore.pyqtSlot
	Property = QtCore.pyqtProperty

def is_pyside():
    if USES_PYSIDE:
        return True
    else:
        return False

if __name__ == "__main__":
	pass