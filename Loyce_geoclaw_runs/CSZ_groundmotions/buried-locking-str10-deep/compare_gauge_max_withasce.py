
import sys, os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

root_dir = os.environ['CHT']

from pylab import *

import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import topotools
from clawpack.visclaw import plottools, geoplot, gaugetools
from clawpack.geoclaw import fgout_tools

event = os.path.split(os.getcwd())[-1]
event = '%s, Varying DEM and/or Grid  Resolution' % event
#asce_region = 1; y1region = 47.3; y2region = 48.6
asce_region = 2; y1region = 46.335; y2region = 47.3
#asce_region = 3; y1region = 45.445; y2region = 46.335
#asce_region = 4; y1region = 44.5; y2region = 45.445
#asce_region = 5; y1region = 43.2; y2region = 44.5
#asce_region = 6; y1region = 42.1; y2region = 43.2
#asce_region = 7; y1region = 41.4; y2region = 42.1
#asce_region = 8; y1region = 40.1; y2region = 41.4

zoom = True
plot_gauge_timeseries = False

#outdir1 = '_output60DEM_15Grid'
outdir1 = None

#outdir2 = '_output15DEM_15Grid_again'
#outdir2 = '_output15DEM_15Grid_large15DEM'
outdir2 = '_output_15sec'

#outdir3 = '_output15GEBDEM_15Grid'
outdir3 = None

#outdir4 = '_output60DEM_30Grid'
outdir4 = None

#outdir5 = '_output15DEM_30Grid'
outdir5 = None

#outdir6 = '_output15DEM_5Grid_large15DEM'
outdir6 = '_output_with5sec'
#outdir6 = None

outdir_list = [outdir1,outdir2,outdir3,outdir4,outdir5,outdir6]
color_list = ['r','b','g','r','m','y']
linestyle_list = ['-','-','-','--','--','--']

fname = root_dir + '/info/asce_values.txt'
d = loadtxt(fname,skiprows=1)
y_asce = d[:,2]
eta_asce = d[:,4]  #*0.3048  # convert from feet to meters

#plotdir = '_plots_with5sec'
#plotdir = '_plots_15Grid_vs_5Grid_large15DEM'
#plotdir = '_plots_15Grid_vs_5Grid_large15DEM'
plotdir = '_plots_5vs15'
os.system('mkdir -p %s' % plotdir)

fname = root_dir + '/topo/topofiles/etopo22_15s_-137_-121_37_55.asc'
etopo = topotools.Topography(fname,3)

plot_eta = True

#gaugenos = range(1,642,20)

def read_gauges(outdir, gaugenos=None):
    if gaugenos is None:
        setgauges = gaugetools.read_setgauges(outdir)
        gaugenos = setgauges.gauge_numbers
    
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
    #print('+++ gaugenos sorted S to N:\n  ',gaugenos)
    gmax = array(gmax)
    xg = array(xg)
    yg = array(yg)
    xg = xg[isort]
    yg = yg[isort]
    gmax = gmax[isort]

    return xg,yg,gmax,gaugenos

    
outdir = outdir2
xg,yg,gmax,gaugenos = read_gauges(outdir)

#============================
# side-by-side plots

if 0:
    figure(201, figsize=(11,8))
    clf()
    ax = axes([0.5,0.1,0.4,0.8])
    ax.contourf(etopo.X, etopo.Y, etopo.Z, [-1e4,0,1e4], 
             colors=[geoplot.blue_green,geoplot.dark_green])


    ax.set_aspect(1./cos(48*pi/180))

    if plot_eta:
        fgframeno = 17
        # Instantiate object for reading fgout frames:
        format = 'binary'  # format of fgout grid output
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
    if zoom:
        ylimits = [44.5,48]
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
     

    lw = 0.8
    for k,outdir in enumerate(outdir_list):
        if outdir is not None:
            print('Reading gauges from ',outdir)
            xg,yg,gmax,gaugenos = read_gauges(outdir)
            base = zeros(len(gmax))
            ampl = base + gmax
            try:
                label = open(outdir + '/label.txt').readline().strip()
            except:
                label = outdir.replace('_','')
            ax2.plot(ampl,yg,color=color_list[k], linestyle=linestyle_list[k],
                      linewidth=lw, label=label)
        

    ax2.set_title('Maximum Amplitude at 100m-depth Gauges\n%s' % event)
    ax2.grid(True)
    ax2.set_xlim(15,0)
    ax2.set_xlabel('meters')
    ax2.set_ylim(ylimits)
    ax2.legend(loc='upper left')

    fname = '%s/map_gauge_comparisons.png' % plotdir
    savefig(fname)
    print('Created %s' % fname)


#==================================
# Rotated to horizontal orientation:

if zoom:
    figno = 203
else:
    figno = 202
figure(figno, figsize=(10,8))
clf()
ax = axes([0.1,0.1,0.8,0.25])
ax.contourf(etopo.Y.T, flipud(etopo.X.T), flipud(etopo.Z.T), [-1e4,0,1e4], 
         colors=[geoplot.blue_green,geoplot.dark_green])


#ax.set_aspect(cos(48*pi/180))

if plot_eta:
    fgframeno = 17
    # Instantiate object for reading fgout frames:
    format = 'binary'  # format of fgout grid output
    fgout_grid = fgout_tools.FGoutGrid(1, outdir, format)      
    fgout1 = fgout_grid.read_frame(fgframeno)
    tfg = fgout1.t / 60.
    eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

    eta_plot = plottools.pcolorcells(fgout1.Y.T,flipud(fgout1.X.T),
                                     flipud(eta.T),
                                     cmap=geoplot.tsunami_colormap)
    clim(-4,4)
    #cb = colorbar(eta_plot, extend='both', shrink=0.5)
    #cb.set_label('meters')
else:
    # plot a few contours of bathymetry:
    CS = ax.contour(etopo.Y.T, flipud(etopo.X.T), flipud(etopo.Z.T),
                    [-1000,-200],
                    linestyles='-', linewidths=0.5, colors='gray')
    #clabel(CS)
    
ax.grid(linewidth=0.3,color='w')
ax.plot(yg,xg,'ro',markersize=1)  # gauge locations
if 0:
    for k in range(3,len(xg),5):
        text(yg[k],xg[k]+0.02,str(gaugenos[k]),fontsize=5,ha='center',va='bottom')

#ax.text(46,-128,'ASCE 100m-depth Gauge Locations',color='r',fontsize=15,
#       ha='center')
ax.text(46.15,-122.2,'Columbia\nRiver',color='w',fontsize=8,
            ha='left',va='bottom')

xlimits = [-122,-130]
ylimits = [40,51]
if zoom:
    ylimits = [44.5,48]
ax.set_xlim(ylimits)
ax.set_ylim(xlimits)
if plot_eta:
    #ax.set_title('Tsunami at %s minutes\nand location of gauges' % tfg)
    ax.text(46,-128,'Tsunami at %s minutes and location of gauges' % tfg,color='r',fontsize=14,
    ha='center')
else:
    ax.set_title('Location of gauges')
ax.set_xlabel('latitude')
ax.set_ylabel('longitude')
    
ax2 = axes([0.1,0.4,0.8,0.5], sharex=ax)
#ax2 = axes([0.75,0.1,0.2,0.8], sharey=ax)

               
lw = 0.8
for k,outdir in enumerate(outdir_list):
    if outdir is not None:
        print('Reading gauges from ',outdir)
        xg,yg,gmax,gaugenos = read_gauges(outdir)
        base = zeros(len(gmax))
        ampl = base + gmax
        try:
            label = open(outdir + '/label.txt').readline().strip()
        except:
            label = outdir.replace('_','')
        ax2.plot(yg,ampl,color=color_list[k], linestyle=linestyle_list[k],
                  linewidth=lw, label=label)


        
ax2.set_title('Maximum Amplitude at Gauges\n%s' % event)
ax2.grid(True)
ax2.set_xlim(ylimits)
#ax2.set_xticks([])
ax2.set_ylabel('meters')
ax2.set_ylim(0,16)

ax2.plot(y_asce,eta_asce,'k-.',lw=lw,markersize=3,label='ASCE 2500-year values')
ax2.legend(loc='upper left', framealpha=1)

#ax2.plot([y1region,y2region],[1,1],'k')
ax2.fill_between([y1region,y2region],[0,0],[16,16],color=[.9,.9,.9])
ax2.text(0.5*(y1region+y2region),0.4,'ASCE\nRegion %i' % asce_region,
         va='bottom',ha='center')

fname = '%s/rotated_map_gauge_comparisons.png' % plotdir
if zoom:
    fname = '%s/rotated_map_gauge_comparisons_zoom.png' % plotdir
savefig(fname)
print('Created %s' % fname)


if plot_gauge_timeseries:
    #=================================
    # compare time series at gauges

    def read_gauge(outdir, gaugeno):
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        t = gauge.t
        eta = gauge.q[3,:]
        xg = gauge.location[0]
        yg = gauge.location[1]
        return t, eta, xg, yg
        
    for gaugeno in gaugenos:
        figure(205,figsize=(12,6))
        clf()
        for k,outdir in enumerate(outdir_list):
            if outdir is None:
                continue
            t,eta,xg,yg = read_gauge(outdir, gaugeno)
            try:
                label = open(outdir + '/label.txt').readline().strip()
            except:
                label = outdir.replace('_','')
            gmax = eta.max()
            label = 'max = %.1fm, ' % gmax + label
            plot(t/60., eta, color=color_list[k], 
                    linestyle=linestyle_list[k], label=label)
        grid(True)
        xlim(0,60)
        xlabel('time (minutes)')
        ylabel('surface elevation (m)')
        legend(framealpha=1)
        title('Gauge %s at (%.4f, %.4f)' % (gaugeno,xg,yg))
        ###
        ### Loyce's note
        ### eta_asce[gaugeno-1] would be the asce value for this gaugeno. 
        ### eta_asce starts with gaugeno = 1 in location eta_asce[0] and
        ### goes sequentially.  It would be good for this ASCE SIFT Region 2 source
        ### to put this value in the title like:
        ### If this gaugeno <= 1620 it can be an ASCE gauge, so include in the title.
        if (gaugeno <= 1620):
            asceval=eta_asce[gaugeno-1]        
            title('Gauge %s at (%.4f, %.4f) had ASCE 2500yr. amp.: %.4f' % (gaugeno,xg,yg,asceval))
        else:
            title('Gauge %s at (%.4f, %.4f)' % (gaugeno,xg,yg))
        ###

        
        fname = '%s/gauge%s_comparison.png' \
                    % (plotdir,str(gaugeno).zfill(5))
        savefig(fname)
        print('Created ',fname)
