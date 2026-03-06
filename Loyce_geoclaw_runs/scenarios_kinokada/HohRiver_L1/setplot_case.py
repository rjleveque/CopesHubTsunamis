"""
setplot function with case argument for looping over many events from
run_geoclaw_make_plots_dtopos.py

This version makes no frame plots, only does post-processing
by calling other script(s).
"""

import os,sys
import matplotlib.pyplot as plt
import clawpack.pyclaw.gauges as gauges
import time, datetime

from pylab import *
from clawpack.geoclaw import topotools

import plot_gaugereport   # Loyces gauge report
import process_fgmax      # Loyces process_fgmax
import make_fgout_Hoh     # Loyces make_fgout_SBIT

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

# MYD directory for this run:
rundir = os.getcwd()

location = 'Hoh' # set to site location if desired to use in plots/filenames

# for a single dtopo run, need to set case dictionary here:
# set here for TACC:
CHTshare = '/work2/04137/rjl/CHTshare/CopesHubTsunamis'
case_single['outdir'] = rundir.replace('/home1/','/scratch/')
case_single['plotdir'] = plotdir.replace('/home1/','/scratch/')
case_single['dtopofile'] = f'{CHTshare}/dtopo/CSZ_Tshirts/CSZ_L1-extended-pmel.tt3'



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
fig_extent[2] = (-127.00,-123.80,45.88,48.0)
var_minmax[2] = (-10.,10.)
topominmax[2] = (-100., 100.)
shrink_colorbar[2] = get_cbshrink(fig_extent[2])

fig_name[3] = 'Hoh 3" to 1" Area'
fig_extent[3] = (-124.5, -124.367, 47.725, 47.767)
var_minmax[3] = (-8.,8.)
topominmax[3] = (-10.,25.)
shrink_colorbar[3] = get_cbshrink(fig_extent[3])

fig_name[4] = 'Hoh 1/3" Area'
fig_extent[4] = (-124.45,-124.3725, 47.725, 47.7525)
var_minmax[4] = (-8.,8.)
topominmax[4] = (-10.,25.)

shrink_colorbar[4] = get_cbshrink(fig_extent[4])

fig_name[5] = 'Hoh Fgmax 1/3" Area'
fig_extent[5] = (-124.435,-124.4025,47.74,47.75)
var_minmax[5] = (-8.,8.)
topominmax[5] = (-10.,25.)
shrink_colorbar[5] = get_cbshrink(fig_extent[5])

fig_name[6] = 'Ruby Beach 1/3" Area'
fig_extent[6] = (-124.417,-124.41, 47.7075, 47.7125)
var_minmax[6] = (-8.,8.)
topominmax[6] = (-10.,25.)
shrink_colorbar[6] = get_cbshrink(fig_extent[6])

fig_name[7] = 'Kalaloch 1/3" Area'
fig_extent[7] = (-124.377,-124.3675, 47.6, 47.6075)
var_minmax[7] = (-8.,8.)
topominmax[7] = (-10.,25.)
shrink_colorbar[7] = get_cbshrink(fig_extent[7])


nfigs=[1,2,3,4,5,6,7]

plt_sites = {}
#This is Gauge 7 at the Chief Klia Wellness Center
plt_sites[1] = ('7', -124.41314815, 47.74101852)

print('\n============= setplot.py ==========================')
print('Start date & time: ', datetime.datetime.now())


def setplot(plotdata, case=case_single):

    """
    Run the code for post-processing.
    This version does not return plotdata for making time frame plots.
    You can also add code to set plotdata and return it if desired.

    case is a dictionary that will be passed in, and it is assumed that
    this dictionary has at least the keys 'outdir', 'plotdir', 'dtopofile'
    with the appropriate values for one particular case (dtopofile).
    """

    from pathlib import Path

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    # these parameters should be set properly when calling this setplot:

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofile = case['dtopofile']

    #assert rundir == case['rundir'], \
    #      f"*** Unexpected: rundir={rundir}\n case['rundir']={case['rundir']}"

    other_figures_dir    = plotdir + '/_other_figures'
    if not os.path.isdir(other_figures_dir):
        os.mkdir(other_figures_dir)
    print('other_figures_dir is: ',other_figures_dir)

    # extract event name:
    if dtopofile is not None:
        event = Path(dtopofile).stem  # drop path and .dtt3 extension
    else:
        event = 'NO_DTOPO'

    run_name = '%s_%s' % (location,event)
    print('In setplot: run_name = ',run_name)

    if 1:
        # call post-processing plot_gaugereport. This prints a .csv and a .txt file:
        plot_gaugereport.report(outdir, plotdir, location, event, dtopofile, run_name)
        print('plot_gaugereport finished')

    if 1:
        # call post-processing process_fgmax for fgno=1
        fgno=1
        process_fgmax.fgreport(fgno,rundir, outdir, plotdir, location, event, dtopofile, run_name)
        print('process_fgmax finished for fgno=1')

        # call post-processing process_fgmax for fgno=2
        fgno=2
        process_fgmax.fgreport(fgno,rundir, outdir, plotdir, location, event, dtopofile, run_name)
        print('process_fgmax finished for fgno=2')

        # call post-processing process_fgmax for fgno=3
        fgno=3
        process_fgmax.fgreport(fgno,rundir, outdir, plotdir, location, event, dtopofile, run_name)
        print('process_fgmax finished for fgno=3')

    if 1:
        # call the fgout animation
        make_fgout_Hoh.make_anim(outdir,plotdir,location,event)
        print('make_fgout_Hoh finished')

    # Make time frames
    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()
    plotdata.format = 'binary'
    plotdata.outdir = outdir

    from clawpack.visclaw import gaugetools
    setgauges = gaugetools.read_setgauges(plotdata.outdir)
    gaugenos = setgauges.gauge_numbers
    print('+++ gaugenos: ', gaugenos)

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
            title('%s\n %s pquake' % (event,timestr))

        def addgauges(current_data,ifig=ifig):
            from clawpack.visclaw import gaugetools
            X1 = fig_extent[ifig][0]; X2 = fig_extent[ifig][1];
            Y1 = fig_extent[ifig][2]; Y2 = fig_extent[ifig][3]

            plt_gnos = []
            for gaugeno in gaugenos:
                G = gauges.GaugeSolution(gaugeno, outdir)
                xplt,yplt = G.location
                if X1 < xplt < X2 and Y1 < yplt < Y2:
                    plt_gnos.append(gaugeno)

            gaugetools.plot_gauge_locations(current_data.plotdata,\
                 gaugenos=plt_gnos,format_string='ko',markersize=0.5,\
                 fontsize=4.,add_labels=True)

        def fixup(current_data,ifig=ifig):
            import pylab
            title_hours(current_data)
            pylab.xticks(fontsize=10)
            pylab.yticks(fontsize=10)
            pylab.ticklabel_format(style='plain',useOffset=False)
            pylab.xticks(rotation=20)
            pylab.gca().set_aspect(1./pylab.cos(47.74 * pylab.pi/180.))
            t = current_data.t

            ### Gauge 7 is the Chief Klia Wellness Center
            if (ifig < 6):
                gaugetools.plot_gauge_locations(current_data.plotdata,gaugenos=[7], \
                    format_string='ko',markersize=1.0,fontsize=6.,add_labels=True)

            if t == 0:
                #addgauges(current_data,ifig=ifig)
                for isite in range(1,len(plt_sites)+1):
                    xplt=plt_sites[isite][1]; yplt=plt_sites[isite][2];
                    site_lbl=plt_sites[isite][0];
                    X1 = fig_extent[ifig][0]; X2 = fig_extent[ifig][1];
                    Y1 = fig_extent[ifig][2]; Y2 = fig_extent[ifig][3]
                    if (ifig < 6):
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

    fname_kmz = 'fgmax_results_Hoh_%s.kmz' % event
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='_other_figures/Hoh_kmlfiles/%s' % fname_kmz)
    fname_kmz = 'fgmax_results_Ruby_%s.kmz' % event
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='_other_figures/Ruby_kmlfiles/%s' % fname_kmz)
    fname_kmz = 'fgmax_results_Kalaloch_%s.kmz' % event
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='_other_figures/Kalaloch_kmlfiles/%s' % fname_kmz)


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

