
import sys, os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *

import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import topotools
from clawpack.visclaw import plottools, geoplot, gaugetools
from clawpack.geoclaw import fgout_tools

location = 'offshore'

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Set CHT enviornment variable to repository top")

etopo_file = os.path.join(CHT,
                          'topo/topofiles/etopo22_15s_-137_-121_37_55.asc')

etopo = topotools.Topography(etopo_file)

plot_eta = True
outdir_instant = None  # not set up for comparisons with instataneous dtopo

def read_gauges(outdir, gaugenos='all'):
    
    if gaugenos == 'all':
        setgauges = gaugetools.read_setgauges(outdir)
        gaugenos = setgauges.gauge_numbers
    
    print('+++ %i gauges in %s' % (len(gaugenos),outdir))
    xg = []
    yg = []
    gmax = []
    for gaugeno in gaugenos:
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        eta = gauge.q[3,:]
        gmax.append(eta.max())
        xg.append(gauge.location[0])
        yg.append(gauge.location[1])

    # sort gauges from south to north
    isort = argsort(yg)
    gaugenos = [gaugenos[i] for i in isort]
    print('+++ gaugenos sorted S to N:\n  ',gaugenos)
    gmax = array(gmax)
    xg = array(xg)
    yg = array(yg)
    xg = xg[isort]
    yg = yg[isort]
    gmax = gmax[isort]

    return xg,yg,gmax

def make_gauge_plot(gaugenos, outdir, plotdir, location, event):
    
    os.system('mkdir -p %s' % plotdir)
    xg,yg,gmax = read_gauges(outdir)

    figure(201, figsize=(11,8))
    clf()
    ax = axes([0.5,0.1,0.4,0.8])
    ax.contourf(etopo.X, etopo.Y, etopo.Z, [-1e4,0,1e4], 
             colors=[geoplot.blue_green,geoplot.dark_green])


    ax.set_aspect(1./cos(48*pi/180))

    if plot_eta:
        fgframeno = 17
        # Instantiate object for reading fgout frames:
        format = 'binary32'  # format of fgout grid output
        fgout_grid = fgout_tools.FGoutGrid(1, outdir, format)      
        fgout1 = fgout_grid.read_frame(fgframeno)
        tfg = fgout1.t / 60.
        eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

        eta_plot = plottools.pcolorcells(fgout1.X,fgout1.Y,eta,
                                         cmap=geoplot.tsunami_colormap)
        clim(-4,4)
        cb = colorbar(eta_plot, extend='both', shrink=0.5)
        cb.set_label('meters')
    else:
        # plot a few contours of bathymetry:
        CS = ax.contour(etopo.X, etopo.Y, etopo.Z, [-1000,-200],
                        linestyles='-', linewidths=0.5, colors='gray')
        #clabel(CS)
        
    ax.grid(linewidth=0.3,color='w')
    ax.plot(xg,yg,'ko',markersize=3)  # gauge locations

    xlimits = [-130,-122]
    ylimits = [40,52]
    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)
    if plot_eta:
        ax.set_title('Tsunami at %s minutes\nand location of gauges' % tfg)
    else:
        ax.set_title('Location of gauges')
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
        
    ax2 = axes([0.1,0.1,0.3,0.8])
    #ax2 = axes([0.75,0.1,0.2,0.8], sharey=ax)

    plot_line = False  # plot horiz line at each gauge amplitude

    if outdir_instant is not None:
        xgi,ygi,gmaxi = read_gauges(outdir_instant)
        base = zeros(len(gmaxi))
        ampl = base + gmaxi
        ax2.plot(ampl,ygi,'r-', label='instantaneous')
        if plot_line:
            for k in range(len(ygi)):
                ax2.plot([base[k],ampl[k]], [ygi[k],ygi[k]], 'r-')    


    base = zeros(len(gmax))
    ampl = base + gmax
    ax2.plot(ampl,yg,'b-', label='time-dependent')
    if plot_line:
        for k in range(len(yg)):
            ax2.plot([base[k],ampl[k]], [yg[k],yg[k]], 'b-')
        
    ax2.set_title('maximum amplitude at gauges\n%s' % event)
    ax2.grid(True)
    ax2.set_xlim(14,0)
    ax2.set_xlabel('meters')
    ax2.set_ylim(ylimits)
    ax2.legend(loc='upper left')


    fname = '%s/map_gauge_amplitudes.png' % plotdir
    savefig(fname)
    print('Created %s' % fname)
    
        
if __name__ == '__main__':
    
    event = os.path.split(os.getcwd())[-1]

    outdir = os.path.abspath('./_output')
    plotdir = os.path.abspath('./_plots')

    gaugenos = 'all'
    make_gauge_plot(gaugenos, outdir,plotdir,location,event)
