#!/usr/bin/env python
import os
import glob

flist = glob.glob("*.ui")
for fname in flist:
    os.system("python pyuic.py %s" % fname)
    print(fname,"->",fname.replace(".ui",".py"))