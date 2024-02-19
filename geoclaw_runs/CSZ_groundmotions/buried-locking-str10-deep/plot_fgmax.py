"""
Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os
import numpy
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import geoplot, plottools
from matplotlib import colors

event = os.path.split(os.getcwd())[-1]
outdir = '_output6'
plotdir = '_plots'

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


clines_zeta = list(numpy.linspace(0,3,7)) + [10]
colors = geoplot.discrete_cmap_1(clines_zeta)
plt.figure(1, figsize=(5,8))
plt.clf()
zeta = numpy.where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
#plt.contourf(fg.X,fg.Y,zeta,clines_zeta,colors=colors)
#plt.colorbar()
pc = plottools.pcolorcells(fg.X, fg.Y, zeta,
                            cmap=cmap_eta, norm=norm_eta)
cb = plt.colorbar(pc, extend='max', shrink=0.7)
cb.set_label('meters')
plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')  # coastline
plt.title("Maximum amplitude\n%s" % event)

if 0:
    # plot arrival time contours and label:
    arrival_t = fg.arrival_time/60.  # arrival time in minutes
    clines_t = numpy.linspace(0,60,10)  # minutes
    clines_t_label = clines_t[2::2]  # which ones to label 
    clines_t_colors = ([.5,.5,.5],)
    con_t = plt.contour(fg.X,fg.Y,arrival_t, clines_t,colors=clines_t_colors) 
    plt.clabel(con_t, clines_t_label)
    plt.title("Maximum amplitude / arrival times\n%s" % event)

# fix axes:
plt.ticklabel_format(style='plain',useOffset=False)
plt.xticks(rotation=20)
plt.gca().set_aspect(1./numpy.cos(fg.Y.mean()*numpy.pi/180.))

if not os.path.isdir(plotdir): 
    os.mkdir(plotdir)
fname = os.path.join(plotdir, "%s_amplitude.png" % event)
plt.savefig(fname)
print("Created ",fname)

