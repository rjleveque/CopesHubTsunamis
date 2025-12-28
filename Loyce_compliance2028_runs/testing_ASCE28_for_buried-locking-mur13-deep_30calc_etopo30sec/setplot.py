"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw import topotools
from clawpack.clawutil.data import ClawData
from clawpack.geoclaw import fgmax_tools
fg = fgmax_tools.FGmaxGrid()

from pylab import *
import clawpack.pyclaw.gauges as gauges
import time
import datetime
import numpy.ma as ma
import matplotlib as mpl

sys.path.insert(0,'.')
import params   # set root_dir, big_topo_dir, loc, event, fgmax_extent
                # only using loc, event, fgmax_extent

# This setplot is for the structure specified in params.loc and thinks it is
# for 'OSVES', so check these are consistent
assert params.loc == 'OSVES', '*** unexpected loc = %s' % params.loc
print('The structure recorded in params.loc was: ',params.loc)

###Note:  fixup below is using Gauge 38 as the closest
###       centered gauge to the OSVES structure

# The event is the earthquake source name used in params.event
print('Using event %s' % params.event)


try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")

print('\n==================================== setplot.py ====================================')
print('Start date & time: ', datetime.datetime.now())

rundir = os.getcwd()
print('rundir = %s' % rundir)

# Set these directories where input data is found:

scenario     = os.path.basename(rundir)
print('scenario is: ',scenario)

plot_dir = params.plot_dir
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)
other_figures_dir    = plot_dir + '/_other_figures'
if not os.path.isdir(other_figures_dir):
    os.mkdir(other_figures_dir)
print('other_figures_dir is: ',other_figures_dir)

fgmax_extent = params.fgmax_extent
print('fgmax_extent was: ',fgmax_extent)
print(' ')


def get_cbshrink(extent):
    X1,X2,Y1,Y2 = extent
    yavg = (Y1+Y2)/2.; cosyavg=cos(yavg*pi/180.)
    LY = 111.319*(Y2-Y1); LX = cosyavg*111.319*(X2-X1); LYdLX = LY/LX; LXdLY = LX/LY
    if LY >= LX: cbshrink = 1.0
    if LY <  LX: cbshrink = LYdLX*7.7/13.
    print('X1,X2,Y1,Y2,LX,LY,LYdLX,LXdLY,cbshrink:', X1,X2,Y1,Y2,LX,LY,LYdLX,LXdLY,cbshrink)
    return cbshrink

fig_name = {}
fig_extent = {}
var_minmax = {}
topominmax ={}
shrink_colorbar = {}


fig_name[1] = 'CompDom'
sec16 = 1./(6*3600.)
cwest = -136.75 - sec16; ceast = -122.75 - sec16;
csouth = 38 - sec16; cnorth = 54 - sec16;
fig_extent[1] = (cwest, ceast, csouth, cnorth)
var_minmax[1] = (-10.,10.)
topominmax[1] = (-100., 100.)
shrink_colorbar[1] = get_cbshrink(fig_extent[1])

fig_name[2] = 'Coast'
fig_extent[2] = (-127.00,-123.80,45.88,47.6)
if 1:  #Want to see better
    var_minmax[2] = (-10.,10.)
topominmax[2] = (-100., 100.)
shrink_colorbar[2] = get_cbshrink(fig_extent[2])

fig_name[3] = 'Ocean Shores'
fig_extent[3] = (-124.3, -124.0, 46.8, 47.1)
if 1: #Want to see better
    var_minmax[3] = (-8.,8.)
topominmax[3] = (-10.,10.)
shrink_colorbar[3] = get_cbshrink(fig_extent[3])

fig_name[4] = 'OSVES'
fig_extent[4] = (-124.1850,-124.1325, 46.9725, 46.99)
if 1:  #Want to see better
    var_minmax[4] = (-8.,8.)
topominmax[4] = (-10.,10.)
shrink_colorbar[4] = get_cbshrink(fig_extent[4])

nfigs=[1,2,3,4]

plt_sites = {}
#plt_sites[1] = ('OSVES', -124.1546, 46.9827)

## This is Option 1, Ocean Shores Elementary School
plt_sites[1] = ('OSVES', -124.154907407, 46.978796296)
#plt_sites[1] = ('CMH', -123.819, 46.188)


#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  
    plotdata.format = 'binary'

    from clawpack.visclaw import gaugetools
    setgauges = gaugetools.read_setgauges(plotdata.outdir)
    gaugenos = setgauges.gauge_numbers
    print('+++ gaugenos: ', gaugenos)

    ### Not using yet
    blue_yellow_red = colormaps.make_colormap({0:[0,0,1], \
                      0.5:[1,.8,.2], 1:[1,0,0]})

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------

    for ifig in nfigs:
        plotfigure = plotdata.new_plotfigure(name=fig_name[ifig], figno=ifig)

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes('pcolor')
        plotaxes.title = fig_name[ifig]

        #plotaxes.scaled = True  #True does not allow aspect ratio

        def timeformat(t):
            # Format the time as hh:mm:ss.ss
            from numpy import mod, sign
            signt = sign(t)
            t = abs(t)
            hours = int(t/3600.) # seconds to integer number of hours
            tmin = mod(t,3600.)  # seconds of remaining time beyond integer number of hours
            min = int(tmin/60.)  # seconds to integer number of minutes
            sec = int(mod(tmin,60.)) # remaining integer sec
            tenth_sec = int(10*(t - int(t)))
            timestr = '%s:%s:%s.%s' % (hours,str(min).zfill(2),str(sec).zfill(2),\
                                       str(tenth_sec).zfill(1))
            if signt < 0:
                timestr = '-' + timestr
            return timestr

        def title_hours(current_data):
            # Title for each frame with time in format hh:mm:ss.ss
            from pylab import title
            t = current_data.t
            timestr = timeformat(t)
            title('%s\n %s pquake' % (scenario,timestr))

        def addgauges(current_data,ifig=ifig):
            from clawpack.visclaw import gaugetools
            X1 = fig_extent[ifig][0]; X2 = fig_extent[ifig][1];
            Y1 = fig_extent[ifig][2]; Y2 = fig_extent[ifig][3]

            plt_gnos = []
            for gaugeno in gaugenos:
                G = gauges.GaugeSolution(gaugeno, '_output')
                xplt,yplt = G.location
                if X1 < xplt < X2 and Y1 < yplt < Y2:
                    plt_gnos.append(gaugeno)

            gaugetools.plot_gauge_locations(current_data.plotdata, \
                 gaugenos=plt_gnos,format_string='ko',markersize=0.5,fontsize=4.,add_labels=True)

        def fixup(current_data,ifig=ifig):
            import pylab
            title_hours(current_data)
            pylab.xticks(fontsize=10)
            pylab.yticks(fontsize=10)
            pylab.ticklabel_format(style='plain',useOffset=False)
            pylab.xticks(rotation=20)

            #Latitude around 47: plt_sites[1] = ('OSVES', -124.154907407, 46.978796296)
            pylab.gca().set_aspect(1./pylab.cos(47.0 * pylab.pi/180.))
            t = current_data.t

            ### Gauge 38 is the closest centered gauge to the OSVES structure
            gaugetools.plot_gauge_locations(current_data.plotdata,gaugenos=[38], \
                format_string='ko',markersize=1.0,fontsize=6.,add_labels=True)

            addgauges(current_data,ifig=ifig)
            if t == 0:
                #addgauges(current_data,ifig=ifig)
                for isite in range(1,len(plt_sites)+1):
                    xplt=plt_sites[isite][1]; yplt=plt_sites[isite][2]; site_lbl=plt_sites[isite][0];
                    X1 = fig_extent[ifig][0]; X2 = fig_extent[ifig][1];
                    Y1 = fig_extent[ifig][2]; Y2 = fig_extent[ifig][3]
                    if X1 < xplt < X2 and Y1 < yplt < Y2:
                        plt.plot(plt_sites[isite][1],plt_sites[isite][2],'k.',markersize=0.5)
                        plt.text(plt_sites[isite][1],plt_sites[isite][2],plt_sites[isite][0]+' ',\
                        fontweight='normal',fontsize=4.,ha='right')
        plotaxes.afteraxes = fixup

        # Water
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        var_min=var_minmax[ifig][0]; var_max=var_minmax[ifig][1]
        plotitem.pcolor_cmin = var_min
        plotitem.pcolor_cmax = var_max
        plotitem.add_colorbar = True
        fig = plt.gcf()
        plotitem.colorbar_shrink = shrink_colorbar[ifig]
        plotitem.plot_var = geoplot.surface_or_depth
        cm_water = geoplot.tsunami_colormap
        plotitem.pcolor_cmap = cm_water
        plotitem.colorbar_label = 'h+B when B<0; h if B>0 (m)'
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]

        # Land
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.land
        Bmin=topominmax[ifig][0]; Bmax=topominmax[ifig][1]
        cm_land = geoplot.land_colors
        plotitem.pcolor_cmap = cm_land
        plotitem.pcolor_cmin = 0
        plotitem.pcolor_cmax = 100
        x1 = fig_extent[ifig][0]; x2 = fig_extent[ifig][1];
        y1 = fig_extent[ifig][2]; y2 = fig_extent[ifig][3]
        plotaxes.xlimits = [x1,x2]
        plotaxes.ylimits = [y1,y2]

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------


    time_scale = 1./3600.
    time_label = 'hours'

    plotfigure = plotdata.new_plotfigure(name='Gauges',figno=300,type='each_gauge')
    plotfigure.kwargs = {'figsize':(8,8)}
    # plotfigure.kwargs = {'figsize':(12,12)}
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(4,1,1)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Quantities of interest'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    # Plot eta, the amount above MHW everywhere
    def eta(current_data):
        t = current_data.t/60.  # Time in minutes
        q = current_data.q
        eta = q[3,:]
        return eta
    plotitem.plot_var = eta
    plotitem.plotstyle = 'b-'

    # Plot h:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(4,1,2)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def h(current_data):
        t = current_data.t/60.  # Time in minutes
        q = current_data.q
        h = q[0,:]
        return h
    plotitem.plot_var = h
    plotitem.plotstyle = 'b-'

    # Plot current speed:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(4,1,3)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def speed(current_data):
        from numpy import where, sqrt
        t = current_data.t/60.  # Time in minutes
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = current_data.q[1,:] / h
        v = current_data.q[2,:] / h
        s = sqrt(u**2 + v**2)
        return s
    plotitem.plot_var = speed
    plotitem.plotstyle = 'b-'

    # Plot momentum flux:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(4,1,4)'
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def mflux(current_data):
        from numpy import where, sqrt
        t = current_data.t/60.  # Time in minutes
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        hu = current_data.q[1,:]
        hv = current_data.q[2,:]
        mf = (hu**2 + hv**2)/h
        return mf
    plotitem.plot_var = mflux
    plotitem.plotstyle = 'b-'

    def aa(current_data):

        from pylab import clf, subplot
        # t = current_data.t  # Time in seconds
        # t = current_data.t/60.  # Time in minutes
        subplot(411)
        ylabel('Surface (eta=h+B), (m)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)

        subplot(412)
        ylabel('h (m)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)

        subplot(413)
        xlabel('')
        ylabel('s, speed (m/s)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)

        subplot(414)
        xlabel('time (hours after earthquake)')
        ylabel('hss, mf (m^3/s^2)')
        title('')
        grid(color='k', linestyle='dotted', linewidth=0.5)

    plotaxes.afteraxes = aa


    #-------------------------------
    # Other Figures for this Site 
    #-------------------------------

    print('other_figures_dir = ',other_figures_dir)
    if os.path.isdir(other_figures_dir):
        files = os.listdir(other_figures_dir)
        print('files: ',files,' END of files')
        if len(files) > 0:
            files.sort()
            print('%i files in this directory' % len(files))
            for filename in files:
                print('\nfilename=',filename)
                path='_other_figures/'+filename
                otherfigure = plotdata.new_otherfigure(name=filename,fname=path)
                print('Added other figure: ',path)
        else:
            print('No files in this directory')
    else:
        print('*** directory not found, will not add to index')

    #-----------------------------------------
    # Figures for fgmax plots
    #-----------------------------------------
    # Note: You need to move fgmax png files into _plots/fgmax_plots after 
    # creating them, e.g., by running the process_fgmax notebook or script.
    # The lines below just create links to these figures from _PlotIndex.html 

    if 0:
        # included by listdir version above:
        otherfigure = plotdata.new_otherfigure(name='max depth',
                        fname='fgmax_plots/h_onshore.png')

        otherfigure = plotdata.new_otherfigure(name='max speed',
                        fname='fgmax_plots/speed.png')
            
    fname_kmz = 'fgmax_results_%s_%s.kmz' % (params.loc,params.event)
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='_other_figures/kmlfiles/%s' % fname_kmz)
        

    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        import os
        from clawpack.visclaw import plot_timing_stats
        try:
            timing_plotdir = plotdata.plotdir + '/timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            units = {'comptime':'hours', 'simtime':'hours', 'cell':'billions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, make_pngs=True,
                                          plotdir=timing_plotdir, units=units)
            os.system('cp %s/timing.* %s' % (plotdata.outdir, timing_plotdir))
        except:
            print('*** Error making timing plots')

    # create a link to this webpage from _PlotIndex.html:
    otherfigure = plotdata.new_otherfigure(name='timing',
                    fname='timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots


    #---------------------------------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                   # print figures
    plotdata.print_format = 'png'               # file format

    # ALL frames and gauges
    plotdata.print_framenos = 'all'             # list of frames to print
    plotdata.print_gaugenos = 'all'             # list of gauges to print
    plotdata.print_fignos   = 'all'             # list of figures to print
    plotdata.html = True                        # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                       # create latex file of plots?
    plotdata.latex_figsperline = 4              # layout of plots 2
    plotdata.latex_framesperline = 4            # layout of plots 1
    plotdata.latex_makepdf = False              # also run pdflatex?
    plotdata.parallel = True                    # Faster

    return plotdata
