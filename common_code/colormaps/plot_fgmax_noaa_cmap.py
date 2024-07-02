"""
Sample code to plot fgmax with NOAA colormap
"""

from pylab import *
import os,sys
from clawpack.geoclaw import fgmax_tools

CHT = os.environ['CHT']   # assuming environment variable set
sys.path.insert(0, os.path.join(CHT,'common_code'))
import noaa_colormaps
cmap_noaa_def, cmap_noaa_max = noaa_colormaps.load()

fg = fgmax_tools.FGmaxGrid()
fg.outdir = '/Users/rjl/git/clawpack/geoclaw/examples/tsunami/chile2010_fgmax-fgout/_output'
data_file = os.path.join(fg.outdir, 'fgmax_grids.data')
fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
fg.read_output()


figure(1)
clf()
plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')
hwet = where(fg.h > 0.1, fg.h, nan)
zeta = where(fg.B>0, hwet, fg.h+fg.B)
zeta = numpy.where(zeta >= 0., zeta, numpy.nan)
pcolormesh(fg.X,fg.Y,zeta,cmap=cmap_noaa_max)
clim(0,0.5)
#clim(-2,10)
colorbar(extend='max')
gca().set_aspect(1/cos(30*pi/180))
title('fgmax with NOAA colormap')

