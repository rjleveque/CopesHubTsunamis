
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from clawpack.geoclaw import topotools
from six.moves import range
import os,sys

from clawpack.visclaw import gaugetools

from clawpack.visclaw import particle_tools
from clawpack.visclaw import legend_tools

import params

cmax = 1.5
cmin = -cmax

cmax_land = 20.


#--------------------------
def setplot(plotdata=None):
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


    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'


    print('Reading all gauges...')
    gauge_solutions = particle_tools.read_gauges(gaugenos='all', 
                                                 outdir=plotdata.outdir)

    gaugenos_lagrangian = [k for k in gauge_solutions.keys() \
                if gauge_solutions[k].gtype=='lagrangian']
    gaugenos_stationary = [k for k in gauge_solutions.keys() \
                if gauge_solutions[k].gtype=='stationary']

    print('+++ gaugenos_lagrangian: ',gaugenos_lagrangian)
    
    def add_particles(current_data):
        t = current_data.t

        # plot recent path:
        t_path_length = 900.   # length of path trailing particle
        kwargs_plot_path = {'linewidth':1, 'color':'k'}
        particle_tools.plot_paths(gauge_solutions, 
                                  t1=t-t_path_length, t2=t, 
                                  gaugenos=gaugenos_lagrangian, 
                                  kwargs_plot=kwargs_plot_path,
                                  extend='both')

        # plot current location:
        kwargs_plot_point = {'marker':'o','markersize':3,'color':'k'}
        particle_tools.plot_particles(gauge_solutions, t, 
                                      gaugenos=gaugenos_lagrangian, 
                                      kwargs_plot=kwargs_plot_point,
                                      extend='both')  

        # plot any stationary gauges:
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos_stationary, format_string='kx', add_labels=False)
        #kwargs={'loc':'upper left'}
        if 0:
            legend_tools.add_legend(['Lagrangian particle','Stationary gauge'],
                    linestyles=['',''], markers=['o','x'],
                    loc='upper right', framealpha=1, fontsize=10)



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


    def fixup(current_data):
        import pylab
        #addgauges(current_data)
        t = current_data.t
        t = t / 60.  # minutes
        pylab.title('Surface at %4.2f minutes' % t, fontsize=10)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, savefig
        gca().set_aspect(1./cos(48*pi/180.))
        title_hours(current_data)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        
    def aa_particles(current_data):
        aa(current_data)
        if len(gaugenos_lagrangian) > 0:
            add_particles(current_data)


    def speed2d(current_data):
        from pylab import sqrt,where,ma
        q = current_data.q
        h = q[0,:,:]
        h = ma.masked_where(h<0.01, h)
        hs = sqrt(q[1,:,:]**2 + q[2,:,:]**2)
        s = hs/h
        return s

    bounds_speed = np.array([1e-6,0.5,1.0,1.5,2,2.5,3,4.5,6])
    cmap_speed = mpl.colors.ListedColormap([[.9,.9,1],[.6,.6,1],\
                     [.3,.3,1],[0,0,1], [1,.8,.8],\
                     [1,.6,.6], [1,.3,.3], [1,0,0]])

    # Set color for value exceeding top of range to purple:
    cmap_speed.set_over(color=[1,0,1])
    norm_speed = mpl.colors.BoundaryNorm(bounds_speed, cmap_speed.N)

    def up_to_level(max_level):
        """
        Make a list like [1,1,1,0] in the case max_level==3.
        Useful if you want to set amr_data_show to only plot on coarse levels
        """
        return max_level*[1] + [0]
                
    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Computational domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    plotaxes.afteraxes = aa
    #plotaxes.afteraxes = add_particles

    ## Limits below never used for AK, CSZ_L1 or SFL
    #plotaxes.xlimits = [-129.16,-122.16]
    #plotaxes.ylimits = [46.0,51.0]


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0,0,0,0]
    plotitem.amr_data_show = up_to_level(4)

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,0]
    plotitem.amr_data_show = up_to_level(4)


    #-----------------------------------------
    # Figure for Outer Coastal area - eta
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name="Outer Coast eta", figno=8)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize': (8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.scaled = False

    plotaxes.xlimits = [-130.6,-123.23]
    plotaxes.ylimits = [46.03,50.97]

    plotaxes.afteraxes = aa
    #plotaxes.afteraxes = aa_particles

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(4)

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(4)

    #-----------------------------------------
    # Figure for Study Region - eta
    #-----------------------------------------
    
    x1,x2,y1,y2 = params.fgmax_extent

    plotfigure = plotdata.new_plotfigure(name="Study Area eta", figno=10)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (8,7)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.scaled = False

    plotaxes.xlimits = [x1,x2]
    plotaxes.ylimits = [y1,y2]

    plotaxes.afteraxes = aa
    #plotaxes.afteraxes = aa_particles

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(5)

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(5)

    #-----------------------------------------
    # Figure for Study Area - speed
    #-----------------------------------------

    plotfigure = plotdata.new_plotfigure(name="Study Area speed", figno=11)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (8,9)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.scaled = False

    plotaxes.xlimits = [x1,x2]
    plotaxes.ylimits = [y1,y2]

    plotaxes.afteraxes = aa
    #plotaxes.afteraxes = aa_particles


    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = speed2d
    plotitem.pcolor_cmap = cmap_speed
    plotitem.kwargs = {'norm': norm_speed}
    plotitem.add_colorbar = True
    plotitem.colorbar_shrink = 0.7
    plotitem.colorbar_kwargs = {'extend':'max'}
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(6)

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax_land
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.patchedges_show = 0
    plotitem.amr_data_show = up_to_level(6)


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    
    time_scale = 1./3600.
    time_label = 'hours'
    
    plotfigure = plotdata.new_plotfigure(name='gauge depth', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def setglimits_depth(current_data):
        from pylab import xlim,ylim,title,argmax,show,array,ylabel
        gaugeno = current_data.gaugeno
        q = current_data.q
        depth = q[0,:]
        t = current_data.t
        g = current_data.plotdata.getgauge(gaugeno)
        level = g.level
        maxlevel = max(level)

        #find first occurrence of the max of levels used by
        #this gauge and set the limits based on that time
        argmax_level = argmax(level)
        xlim(time_scale*array(t[argmax_level],t[-1]))
        ylabel('meters')
        min_depth = depth[argmax_level:].min()
        max_depth = depth[argmax_level:].max()
        ylim(min_depth-0.5, max_depth+0.5)
        title('Gauge %i : Flow Depth (h)\n' % gaugeno + \
              'max(h) = %7.3f,    max(level) = %i' %(max_depth,maxlevel))    
        #show()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label

    # Plot depth as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    ## Set the limits and the title in the function below
    plotaxes.afteraxes = setglimits_depth

    plotfigure = plotdata.new_plotfigure(name='gauge surface eta', figno=301, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def setglimits_eta(current_data):
        from pylab import xlim,ylim,title,argmax,show,array,ylabel
        gaugeno = current_data.gaugeno
        q = current_data.q
        eta = q[3,:]
        t = current_data.t
        g = current_data.plotdata.getgauge(gaugeno)
        level = g.level
        maxlevel = max(level)

        #find first occurrence of the max of levels used by
        #this gauge and set the limits based on that time
        argmax_level = argmax(level) #first occurrence of it
        xlim(time_scale*array(t[argmax_level],t[-1]))
        ylabel('meters')
        min_eta = eta[argmax_level:].min()
        max_eta = eta[argmax_level:].max()
        ylim(min_eta-0.5,max_eta+0.5)
        title('Gauge %i : Surface Elevation (eta)\n' % gaugeno + \
              'max(eta) = %7.3f,    max(level) = %i' %(max_eta,maxlevel))
        #show()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label
    
    # Plot surface (eta) as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    ## Set the limits and the title in the function below
    plotaxes.afteraxes = setglimits_eta

    plotfigure = plotdata.new_plotfigure(name='speed', figno=302, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    def speed(current_data):
        from numpy import sqrt, maximum, where, diff, hstack, cos, pi
        q = current_data.q
        t = current_data.t
        gaugeno = current_data.gaugeno
        if gaugeno in gaugenos_stationary:
            h   = q[0,:]
            hu  = q[1,:]
            hv  = q[2,:]
            s   = sqrt(hu**2 + hv**2) / maximum(h,0.001)
            s   = where(h > 0.001, s, 0.0)
            
        elif gaugeno in gaugenos_lagrangian:
            print('Computing speed for particle by differencing x,y')
            x = q[1,:]
            y = q[2,:]
            dt = diff(t)
            # difference x,y and convert from degrees to meters:
            u = diff(x)/dt * 111e3 * cos(48*pi/180.)
            v = diff(y)/dt * 111e3
            s = sqrt(u**2 + v**2)
            s = hstack((s[0],s))  # so it's the same length as t
            
        else:
            raise ValueError('Unexpected gaugeno')

        return s

    def setglimits_speed(current_data):
        from pylab import xlim,ylim,title,argmax,show,array,ylabel
        gaugeno = current_data.gaugeno
        s = speed(current_data)
        t = current_data.t
        g = current_data.plotdata.getgauge(gaugeno)
        level = g.level
        maxlevel = max(level)

        #find first occurrence of the max of levels used by
        #this gauge and set the limits based on that time
        argmax_level = argmax(level) #first occurrence of it
        xlim(time_scale*array(t[argmax_level],t[-1]))
        ylabel('meters/sec')
        min_speed = s[argmax_level:].min()
        max_speed = s[argmax_level:].max()
        ylim(min_speed-0.5,max_speed+0.5)
        title('Gauge %i : Speed (s)\n' % gaugeno + \
              'max(s) = %7.3f,    max(level) = %i' %(max_speed,maxlevel))
        #show()

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.time_scale = time_scale
    plotaxes.time_label = time_label

    # Plot speed (s) as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = speed
    plotitem.plotstyle = 'b-'

    ## Set the limits and the title in the function below
    plotaxes.afteraxes = setglimits_speed



    #-----------------------------------------
    # Figures for fgmax plots
    #-----------------------------------------
    # Note: You need to move fgmax png files into _plots/fgmax_plots after 
    # creating them, e.g., by running the process_fgmax notebook or script.
    # The lines below just create links to these figures from _PlotIndex.html 

    otherfigure = plotdata.new_otherfigure(name='max depth',
                    fname='fgmax_plots/h_onshore.png')

    otherfigure = plotdata.new_otherfigure(name='max speed',
                    fname='fgmax_plots/speed.png')
        
    fname_kmz = 'fgmax_results_%s_%s.kmz' % (params.loc,params.event)
    otherfigure = plotdata.new_otherfigure(name=fname_kmz,
                    fname='fgmax_kml/%s' % fname_kmz)
        
    # add additional lines for any other figures you want added to the index.            
                    

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


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    #plotdata.print_framenos = list(range(21,29))     # list of frames to print
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'             # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata

