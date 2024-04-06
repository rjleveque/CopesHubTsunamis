
import sys, os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *

import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import topotools
from clawpack.visclaw import geoplot

event = os.path.split(os.getcwd())[-1]
root_dir = os.environ['CHT']

outdir = '_output'
plotdir = '_plots'
os.system('mkdir -p %s' % plotdir)

topodir = root_dir + '/topo/topofiles'
etopo = topotools.Topography(topodir + '/etopo22_1min_-163_-122_38_63.asc',3)
#etopo = topotools.Topography('/Users/rjl/topo/topofiles/etopo1_-163_-122_38_63.asc',3)


gaugenos = range(1,642,20)
xg = []
yg = []
gmax = []
for gaugeno in gaugenos:
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    eta = gauge.q[3,:]
    gmax.append(eta.max())
    xg.append(gauge.location[0])
    yg.append(gauge.location[1])

gmax = array(gmax)
xg = array(xg)
yg = array(yg)

figure(200, figsize=(12,7))
clf()

plot(yg, gmax, 'k-o')
xlabel('latitude')
ylabel('offshore wave amplitude (m)')
grid(True)
title('Gauge amplitudes\n%s' % event)
fname = '%s/gauge_amplitudes.png' % plotdir
savefig(fname)
print('Created %s' % fname)


figure(201, figsize=(10,8))
clf()
ax = axes([0.6,0.1,0.3,0.8])
ax.contourf(etopo.X, etopo.Y, etopo.Z, [-1e4,0,1e4], 
         colors=[geoplot.blue_green,geoplot.dark_green])

CS = ax.contour(etopo.X, etopo.Y, etopo.Z, [-1000,-200],
                linestyles='-', linewidths=0.5, colors='gray')
#clabel(CS)

ax.set_aspect(1./cos(48*pi/180))
         

ax.plot(xg,yg,'ko',markersize=3)  # gauge locations

xlimits = [-130,-122]
ylimits = [40,52]
ax.set_xlim(xlimits)
ax.set_ylim(ylimits)
ax.set_title('Location of gauges')
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
    
ax2 = axes([0.1,0.1,0.3,0.8])
#ax2 = axes([0.75,0.1,0.2,0.8], sharey=ax)
base = zeros(len(gmax))
ampl = base + gmax
ax2.fill_between(base,ampl,color=[.6,.6,1])


ax2.plot(ampl,yg,'b-')
for k in range(len(yg)):
    ax2.plot([base[k],ampl[k]], [yg[k],yg[k]], 'b-')
        

ax2.set_title('maximum amplitude at gauges\n%s' % event)
ax2.grid(True)
ax2.set_xlim(15,0)
ax2.set_xlabel('meters')
ax2.set_ylim(ylimits)


fname = '%s/map_gauge_amplitudes.png' % plotdir
savefig(fname)
print('Created %s' % fname)
