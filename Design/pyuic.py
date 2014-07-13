#!/usr/bin/env python

import os
import sys
if len(sys.argv) != 2:
    raise TypeError("usage: %s sample.ui" % sys.argv[0])

infile = sys.argv[1]
outfile = infile.replace(".ui", ".py")
os.system("pyuic4-3.3 %s | sed -e 's/PyQt4/PySide/g' > %s" % (infile, outfile))    
os.system('gsed -i -e "s/QtGui.QMenuBar(MainWindow)/QtGui.QMenuBar() #QtGui.QStatusBar(MainWindow) /g" %s' % outfile)
os.system('gsed -i -e "s/QtGui.QMenuBar(MainWindow)/QtGui.QMenuBar() #QtGui.QStatusBar(MainWindow) /g" %s' % outfile)
os.system('gsed -i -e "s/from mplwidget import MplWidget/from .mplwidget import MplWidget/g" %s' % outfile)
