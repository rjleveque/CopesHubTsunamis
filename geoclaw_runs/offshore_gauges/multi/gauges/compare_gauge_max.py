
import sys, os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *

import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import topotools
from clawpack.visclaw import plottools, geoplot, gaugetools
from clawpack.geoclaw import fgout_tools

#event = os.path.split(os.getcwd())[-1]
event = 'Comparison of different events'
location = 'offshore'

zoom = False
plot_gauge_timeseries = False
add_asce_2500 = False
all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


if 1:

    models = all_models[3:]
    subtitle = 'locking models'
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

color_list = 3*['r','b','g','r','b','g','r','b','g']
linestyle_list = 3*['-','-','-','--','--','--','-.','-.','-.']

d = loadtxt('/Users/rjl/git/CopesHubTsunamis/info/asce_values.txt',skiprows=1)
y_asce = d[:,2]
eta_asce = d[:,4]  #*0.3048  # convert from feet to meters

plotdir = '.'

etopo = topotools.Topography('/Users/rjl/topo/topofiles/etopo1_-163_-122_38_63.asc',3)

plot_eta = False

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

    
outdir = '/Users/rjl/git/CopesHubTsunamis/geoclaw_runs/offshore_gauges/multi/scratch/geoclaw_outputs/_output_buried-locking-str10-deep'
xg,yg,gmax,gaugenos = read_gauges(outdir)

#fname = '%s_%s_gauges_max.txt' % (location,events[0])
#yg,gmax = loadtxt(fname, unpack=True)

#==================================
# Rotated to horizontal orientation:

if zoom:
    figno = 203
else:
    figno = 204
figure(figno, figsize=(13,8))
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
    cb = colorbar(eta_plot, extend='both', shrink=0.5)
    cb.set_label('meters')
else:
    # plot a few contours of bathymetry:
    CS = ax.contour(etopo.Y.T, flipud(etopo.X.T), flipud(etopo.Z.T),
                    [-1000,-200],
                    linestyles='-', linewidths=0.5, colors='gray')
    #clabel(CS)
    
ax.grid(linewidth=0.3,color='w')
ax.plot(yg,xg,'ro',markersize=1)  # gauge locations
if zoom:
    for k in range(3,len(xg),5):
        text(yg[k],xg[k]+0.02,str(gaugenos[k]),fontsize=5,ha='center',va='bottom')

ax.text(44,-128,'ASCE 100m-depth Gauge Locations',color='r',fontsize=15)
ax.text(46.15,-122.2,'Columbia\nRiver',color='w',fontsize=8,
            ha='left',va='bottom')

xlimits = [-122,-130]
ylimits = [40,51]
if zoom:
    ylimits = [45.5,48]
ax.set_xlim(ylimits)
ax.set_ylim(xlimits)

ax.set_xlabel('latitude')
ax.set_ylabel('longitude')
    
ax2 = axes([0.1,0.4,0.8,0.5], sharex=ax)
#ax2 = axes([0.75,0.1,0.2,0.8], sharey=ax)

               
lw = 0.8
for k,event in enumerate(events):
    fname = '%s_%s_gauges_max.txt' % (location,event)
    try:
        yg,gmax = loadtxt(fname, unpack=True)
        label = event
        ax2.plot(yg,gmax,color=color_list[k], linestyle=linestyle_list[k],
                  linewidth=lw, label=label)
    except:
        print('Failed to load ',fname)

        
ax2.set_title('Maximum Amplitude at Gauges\n%s' % subtitle)
ax2.grid(True)
ax2.set_xlim(ylimits)
#ax2.set_xticks([])
ax2.set_ylabel('meters')
ax2.set_ylim(0,20)

if add_asce_2500:
    ax2.plot(y_asce,eta_asce,'k-.',lw=lw,markersize=3,
             label='ASCE 2500-year values')
             
ax2.legend(loc='upper right', framealpha=1)

fname = '%s/%s_gauge_comparisons_%s.png' \
        % (plotdir,location,subtitle.replace(' ','-'))
if zoom:
    fname = '%s/%s_gauge_comparisons_%s_zoom.png' \
            % (plotdir,location,subtitle.replace(' ','-'))
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
        
    gaugenos = range(5,651,1)

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
###
###     asceval=eta_asce[gaugeno-1]        
###     title('Gauge %s at (%.4f, %.4f) had ASCE 2500yr. amp.: %.4f' % (gaugeno,xg,yg,asceval))
###
        
        fname = '%s/gauge%s_comparison.png' \
                    % (plotdir,str(gaugeno).zfill(5))
        savefig(fname)
        print('Created ',fname)
