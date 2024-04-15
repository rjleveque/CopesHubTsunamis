
"""
Read in one static displacement file to get the lon lat points used
for all, and also for time-dependent ruptures.

Save as lonlat_points.txt.
This file is needed in make_dtopos.py
"""

from pylab import *
import os

event = 'buried-locking-mur13-shallow'
datadir = 'vertical_displacements'
fname_orig = 'vert_displacements_all_xgrid_' + event
path_orig = os.path.join(datadir, fname_orig)
lon,lat,zdisp = loadtxt(path_orig, skiprows=1,usecols=[1,2,3],unpack=True)
points = vstack((lon,lat)).T
fname = 'lonlat_points.txt'
savetxt(fname,points,fmt=' %.8f')
print('Created ',fname)
