"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

NOTE: Uses some new plot attributes that will be introduced in v5.9.1.

"""
import os,sys
import numpy as np
from clawpack.clawutil.util import fullpath_import

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

plot_gauges_site = fullpath_import(f'{CHT}/common_code/plot_gauges_site.py')

# Change these if copied to a different Run directory:

location = 'WestportBridges'




#--------------------------
def setplot(plotdata=None, case={}):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()
    plotdata.format = 'binary'

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofiles = case['dtopofiles']

    dtopofile = dtopofiles[0][-1]
    event = os.path.splitext(os.path.split(dtopofile)[-1])[0]

    plotdata.outdir = outdir
    plotdata.plotdir = plotdir

    run_name = '%s_%s' % (location,event)
    print('In setplot: run_name = ',run_name)
    print('outdir = %s' % outdir)
    print('plotdir = %s' % plotdir)

    if 1:
        gauges_plotdir = plotdir + '/gauges'
        os.system('mkdir -p %s' % gauges_plotdir)
        plot_gauges_site.make_all_plots_and_report(outdir, gauges_plotdir,
                         location=location, event=event,
                         gaugenos='all', sea_level=0.)

    if 0:
        # NEED TO FIX
        try:
            #fgmax_plotdir = plotdir + '/fgmax'
            fgmax_plotdir = plotdir  # fgmax1 or fgmax2 now added in function
            os.system('mkdir -p %s' % fgmax_plotdir)
            process_fgmax.make_all_fgmax_plots(outdir, fgmax_plotdir,
                                           location=location, event=event)
        except:
            print('*** problem with process_fgmax in %s, skipping' \
                  % run_name)
    if 0:
        make_fgout_animation.make_anim(outdir, plotdir, location, event)

    if 1:
        # to suppress frame plots, quit after plotting gauges, fgmax, fgout:
        plotdata = None
        return plotdata


    # for frame plots:

    if 1:
        from clawpack.visclaw import gaugetools
        setgauges = gaugetools.read_setgauges(plotdata.outdir)
        gaugenos = setgauges.gauge_numbers
        print('+++ gaugenos: ', gaugenos)


    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata,
             gaugenos=gaugenos, format_string='ko', add_labels=True,
             fontsize=8, markersize=3)


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Ocean Surface', figno=0)
    #plotfigure.show = False
    # new options in v5.10.0:
    plotfigure.figsize = (7,7)
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')

    # new option: including h:m:s in title converts to hours:minutes:seconds
    plotaxes.title = 'Surface at time h:m:s after quake'

    # new options:
    plotaxes.aspect_latitude = -20.  # set aspect ratio based on this latitude
    plotaxes.title_fontsize = 15
    plotaxes.xticks_kwargs = {'rotation':20, 'fontsize':12}
    plotaxes.yticks_fontsize = 12
    plotaxes.xlabel = 'Longitude'
    plotaxes.xlabel_fontsize = 12
    plotaxes.ylabel = 'Latitude'
    plotaxes.ylabel_fontsize = 12

    plotaxes.xlimits = [-130,-123]
    plotaxes.ylimits = [38,51]

    #plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -2
    plotitem.pcolor_cmax = 2
    plotitem.add_colorbar = True
    plotitem.colorbar_extend = 'both'
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_label = 'meters'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 500.0
    plotitem.add_colorbar = False


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Westport', figno=1)
    # new options in v5.10.0:
    plotfigure.figsize = (8,6)
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')

    # new option: including h:m:s in title converts to hours:minutes:seconds
    plotaxes.title = 'Surface at time h:m:s after quake'

    # new options:
    plotaxes.aspect_latitude = -20.  # set aspect ratio based on this latitude
    plotaxes.title_fontsize = 15
    plotaxes.xticks_kwargs = {'rotation':20, 'fontsize':10}
    plotaxes.yticks_fontsize = 12
    plotaxes.xlabel = 'Longitude'
    plotaxes.xlabel_fontsize = 12
    plotaxes.ylabel = 'Latitude'
    plotaxes.ylabel_fontsize = 12
    plotaxes.useOffset = False

    plotaxes.xlimits = [-124.175, -123.94]
    plotaxes.ylimits = [46.75, 46.93]

    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -5
    plotitem.pcolor_cmax = 5
    plotitem.add_colorbar = True
    plotitem.colorbar_extend = 'both'
    plotitem.colorbar_shrink = 0.5
    plotitem.colorbar_label = 'meters'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0,0,0,0,0,0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 10.0
    plotitem.add_colorbar = False

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------


    time_scale = 1./60.
    time_label = 'minutes'

    plotfigure = plotdata.new_plotfigure(name='Gauges',figno=300,type='each_gauge')
    plotfigure.figsize = (10,5)
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Water depth h'
    plotaxes.grid = True
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'



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
