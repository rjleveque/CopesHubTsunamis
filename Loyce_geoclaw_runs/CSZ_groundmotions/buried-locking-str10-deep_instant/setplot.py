
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os
import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools
from clawpack.amrclaw import region_tools
from clawpack.visclaw import plottools

event = os.path.split(os.getcwd())[-1]

root_dir = os.environ['CHT']
RRdir = root_dir + '/topo/regions/'

if 0:
    image = plt.imread('GE_PA2.png')

    def background_image(current_data):
        from pylab import imshow
        extent = [-123.5,-123.33,48.10,48.17]
        imshow(image,extent=extent)

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

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('%s after earthquake' % timestr)

    def surface_or_depth_feet(current_data):
        sdf = geoplot.surface_or_depth(current_data)
        sdf *= 1./0.3048  # convert to feet
        return sdf
    

    def plot_RR(current_data):
        from pylab import plot
        RRs = ['RuledRectangle_Coast_40_46',
               'RuledRectangle_Coast_40_46b',
               'RuledRectangle_Coast_46_51',
               'RuledRectangle_Coast_46_51b']
        for RR in RRs:
                RRfile = os.path.abspath(RRdir + '/%s.data' % RR)
                RR = region_tools.RuledRectangle()
                RR.read(RRfile)
                xr,yr = RR.vertices()
                plot(xr,yr,'k')
        # also add a rectangle:
        plottools.plotbox([-125.9,-124.1, 46.587, 47.227])

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Domain', figno=0)
    #plotfigure.show = False
    plotfigure.figsize = (6,9)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    #plotaxes.xlimits = [-127,-122]
    #plotaxes.ylimits = [45,49.5]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)
        ticklabel_format(useOffset=False)
        gca().set_aspect(1./cos(48*pi/180.))
        title_hours(current_data)
        plot_RR(current_data)
    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -4. 
    plotitem.pcolor_cmax = 4. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_extend = 'both'
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='CSZ', figno=1)
    #plotfigure.show = False
    plotfigure.figsize = (8,8)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = False
    plotaxes.xlimits = [-130,-122.5]
    plotaxes.ylimits = [38.5,50.5]

    def fixup(current_data):
        from pylab import title, ticklabel_format, gca, cos, pi
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        title('Surface at %4.2f hours' % t, fontsize=10)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)
        ticklabel_format(useOffset=False)
        gca().set_aspect(1./cos(48*pi/180.))
        title_hours(current_data)
        plot_RR(current_data)

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -4. 
    plotitem.pcolor_cmax = 4. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_extend = 'both'
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land1_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 2000.0
    plotitem.add_colorbar = False
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel
        t = current_data.t 
        gaugeno = current_data.gaugeno


        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('time (hours)')

    plotaxes.afteraxes = add_zeroline

    #-----------------------------------------
    # Figures for fgmax plots
    #-----------------------------------------
    # Note: You need to move fgmax png files into the place indicated below after
    # creating them. This file will be in the plots directory indicated in the Makefile
    # The lines below just create links to these figures from _PlotIndex.html

    if 1:
        # included by listdir version above:
        otherfigure = plotdata.new_otherfigure(name='max amplitude',
                        fname='%s_amplitude.png' %event)

    #fname_kmz = 'fgmax_results_%s_%s.kmz' % (params.loc,params.event)
    #otherfigure = plotdata.new_otherfigure(name=fname_kmz,
    #                fname='fgmax_kml/%s' % fname_kmz)


    #-----------------------------------------



    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True

    return plotdata

