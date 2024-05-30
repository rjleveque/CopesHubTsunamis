
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np
import matplotlib.pyplot as plt

from clawpack.geoclaw import topotools

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
    

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Domain', figno=0)
    #plotfigure.show = False
    plotfigure.figsize = (6,9)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface at h:m:s'
    plotaxes.useOffset = False
    plotaxes.aspect_latitude = 48.
    #plotaxes.xlimits = [-130,-122]
    #plotaxes.ylimits = [40,51]

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
    plotitem.patchedges_show = 1

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
    plotfigure.show = False
    plotfigure.figsize = (8,8)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface at h:m:s'
    plotaxes.scaled = False
    plotaxes.xlimits = [-130,-122.5]
    plotaxes.ylimits = [38.5,50.5]
    plotaxes.useOffset = False
    plotaxes.aspect_latitude = 48.


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
    # Zoom for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Port Angeles Surface', figno=2)
    #plotfigure.show = False
    plotfigure.figsize = (8,7)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface at h:m:s\nPort Angeles'
    plotaxes.scaled = False
    plotaxes.xlimits = [-123.5,-123.36]
    plotaxes.ylimits = [48.1,48.2]
    plotaxes.useOffset = False
    plotaxes.aspect_latitude = 48.
    plotaxes.xticks_kwargs = {'rotation':20}


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    #plotitem.plot_var = surface_or_depth_feet
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -2. 
    plotitem.pcolor_cmax = 2. 
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
    # Zoom for speed
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Port Angeles Speed', figno=3)
    #plotfigure.show = False
    plotfigure.figsize = (8,7)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Speed at h:m:s\nPort Angeles'
    plotaxes.xlimits = [-123.5,-123.36]
    plotaxes.ylimits = [48.1,48.2]
    plotaxes.useOffset = False
    plotaxes.aspect_latitude = 48.
    plotaxes.xticks_kwargs = {'rotation':20}


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.speed
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 3. 
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_extend = 'max'
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

