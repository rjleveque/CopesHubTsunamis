"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.


"""
import os,sys
import numpy as np

# Change these if copied to a different Run directory:


if 0:
    rundir = os.getcwd()
    outdir = os.path.join(rundir,'_output')  # change if necessary
    plotdir = os.path.join(rundir,'_plots')  # change if necessary
    #outdir = '.'
    #plotdir = '.'


    print('rundir = %s' % rundir)
    print('outdir = %s' % outdir)
    print('plotdir = %s' % plotdir)



try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW environment variable")



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


    if 1:
        from clawpack.visclaw import gaugetools
        setgauges = gaugetools.read_setgauges(plotdata.outdir)
        gaugenos = setgauges.gauge_numbers
        print('+++ gaugenos: ', gaugenos)


    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata,
             gaugenos=[1001], format_string='ko', add_labels=True,
             fontsize=8, markersize=3)


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name='Ocean Surface', figno=1)
    #plotfigure.show = False
    # new options in v5.10.0:
    plotfigure.figsize = (7,7)
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')

    # new option: including h:m:s in title converts to hours:minutes:seconds
    plotaxes.title = 'Surface at time h:m:s after quake'

    # new options:
    plotaxes.aspect_latitude = 44.  # set aspect ratio based on this latitude
    plotaxes.title_fontsize = 15
    plotaxes.xticks_kwargs = {'rotation':20, 'fontsize':12}
    plotaxes.yticks_fontsize = 12
    plotaxes.xlabel = 'Longitude'
    plotaxes.xlabel_fontsize = 12
    plotaxes.ylabel = 'Latitude'
    plotaxes.ylabel_fontsize = 12

    #plotaxes.xlimits = [-130,-123]
    #plotaxes.ylimits = [38,51]

    plotaxes.afteraxes = addgauges

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
    plotitem.amr_patchedges_show = [0,0,1,1,0,0]
    plotitem.amr_patchedges_color = 6*['k']
    plotitem.amr_patchedges_color[2] = 'g'
    plotitem.amr_patchedges_color[3] = 'r'

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

    plotfigure = plotdata.new_plotfigure(name='Newport-3sec', figno=2)
    # new options in v5.10.0:
    plotfigure.figsize = (8,6)
    plotfigure.facecolor = 'w'

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')

    # new option: including h:m:s in title converts to hours:minutes:seconds
    plotaxes.title = 'Surface at time h:m:s after quake'

    # new options:
    plotaxes.aspect_latitude = 44.  # set aspect ratio based on this latitude
    plotaxes.title_fontsize = 15
    plotaxes.xticks_kwargs = {'rotation':20, 'fontsize':10}
    plotaxes.yticks_fontsize = 12
    plotaxes.xlabel = 'Longitude'
    plotaxes.xlabel_fontsize = 12
    plotaxes.ylabel = 'Latitude'
    plotaxes.ylabel_fontsize = 12
    plotaxes.useOffset = False

    plotaxes.xlimits = [-124.24,-123.9150]
    plotaxes.ylimits = [44.535,44.715]

    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -10
    plotitem.pcolor_cmax = 10
    plotitem.add_colorbar = True
    plotitem.colorbar_extend = 'both'
    plotitem.colorbar_shrink = 0.5
    plotitem.colorbar_label = 'meters'
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0,0,1,1]
    plotitem.amr_patchedges_color = 6*['k']
    plotitem.amr_patchedges_color[5] = 'r'


    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 10.0
    plotitem.add_colorbar = False
    plotitem.amr_patchedges_show = [0] #[0,0,0,0,0,1,1,1]

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
