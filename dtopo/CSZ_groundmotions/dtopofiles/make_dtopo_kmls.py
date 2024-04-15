
from pylab import *
import os,sys,glob

CHT = os.environ['CHT']
sys.path.insert(0, os.path.join(CHT, 'common_code'))
from dtopo_contours_kml import make_kml

dtopofiles = glob.glob('*_instant.dtt3')
for dtopofile in dtopofiles:
    make_kml(dtopofile)
