"""
Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os,sys
import numpy
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import geoplot, plottools
from matplotlib import colors

CHT = os.environ['CHT']   # assuming environment variable set
sys.path.insert(0, os.path.join(CHT,'common_code'))
import noaa_colormaps
cmap_noaa_def, cmap_noaa_max = noaa_colormaps.load()

#event = os.path.split(os.getcwd())[-1]
event = 'buried-random-mur13-deep'
outdir = '%s_swe_shallow/_output_30sec' % event
#plotdir = '_plots'
plotdir = '_plots_agu24'
fname_png = '%s/amplitude_%s_swe_30sec.png' % (plotdir,event)


dtopofile = os.path.join(CHT,'dtopo/CSZ_groundmotions/dtopofiles', 
                        '%s.dtt3' % event)


try:
    title_event = open(outdir + '/label.txt').readline().strip()
except:
    title_event = outdir.replace('_','-')

if 0:
    bounds_eta = [0.1] + list(numpy.linspace(0.5,12,6))

    cmap_eta = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],                 [1,.7,.7], [1,.4,.4], [1,0,0]])

    # Set color for value exceeding top of range to purple:
    cmap_eta.set_over(color=[1,0,1])

    # Set color for land points without inundation to light green:
    cmap_eta.set_under(color=[1,1,1])

    norm_eta = colors.BoundaryNorm(bounds_eta, cmap_eta.N)


os.system('mkdir -p %s' % plotdir)

fg = fgmax_tools.FGmaxGrid()
fg.outdir = outdir
data_file = os.path.join(outdir, 'fgmax_grids.data')
fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
fg.read_output()

fg.interp_dz(dtopofile, dtopo_type=3)
Bfinal = fg.B + fg.dz  # since fg.B captured at t=0 when dz=0 and hence is B0

plt.figure(1, figsize=(5,8))
plt.clf()
#zeta = numpy.where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
hwet = numpy.where(fg.h > 0.1, fg.h, numpy.nan)
#zeta = numpy.where(fg.B>0, hwet, fg.h+fg.B)
zeta = numpy.where(fg.B>0, hwet, fg.h+Bfinal)
zeta = numpy.where(zeta >= 0., zeta, numpy.nan)
pc = plottools.pcolorcells(fg.X, fg.Y, zeta,cmap=cmap_noaa_max)
plt.clim(-2,10)
cb = plt.colorbar(pc, extend='max', shrink=0.7)
cb.set_label('meters')
plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k',linewidths=0.5)  # coastline
plt.title("Maximum amplitude\n%s\n%s" % (event,title_event))
plt.ylim(40,51)

if 0:
    # plot arrival time contours and label:
    arrival_t = fg.arrival_time/60.  # arrival time in minutes
    clines_t = numpy.linspace(0,60,10)  # minutes
    clines_t_label = clines_t[2::2]  # which ones to label 
    clines_t_colors = ([.5,.5,.5],)
    con_t = plt.contour(fg.X,fg.Y,arrival_t, clines_t,colors=clines_t_colors) 
    plt.clabel(con_t, clines_t_label)
    plt.title("Maximum amplitude / arrival times\n%s" % title_event)

# fix axes:
plt.ticklabel_format(style='plain',useOffset=False)
plt.xticks(rotation=20)
plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))

if not os.path.isdir(plotdir): 
    os.mkdir(plotdir)
plt.savefig(fname_png)
print("Created ",fname_png)
