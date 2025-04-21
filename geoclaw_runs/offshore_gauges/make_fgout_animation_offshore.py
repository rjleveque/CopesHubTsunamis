"""
Make an mp4 animation of fgout grid results. 
This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.

"""

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from datetime import timedelta

from clawpack.geoclaw import fgout_tools
    

location = 'offshore'
fgno = 1  # which fgout grid

format = 'binary'  # format of fgout grid output

#fgframes = range(1,241)  # frames of fgout solution to use in animation
#fgframes = range(1,81)  # 40 minutes
fgframes = range(1,82,2)

def make_anim(outdir, plotdir, location, event):
    figsize = (5,8)

    # Instantiate object for reading fgout frames:
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 
    fgout_grid.read_fgout_grids_data()


    # Plot one frame of fgout data and define the Artists that will need to
    # be updated in subsequent frames:

    fgout1 = fgout_grid.read_frame(fgframes[0])

    #plot_extent = fgout1.extent_edges
    plot_extent = [-130,-122.5,38.5,50.5]

    ylat = fgout1.Y.mean()  # for aspect ratio of plots

    fig,ax = subplots(figsize=figsize)

    ax.set_xlim(plot_extent[:2])
    ax.set_ylim(plot_extent[2:])

    plottools.pcolorcells(fgout1.X,fgout1.Y,fgout1.B, cmap=geoplot.land1_colormap)
    clim(0,2000)

    eta = ma.masked_where(fgout1.h<0.001, fgout1.eta)

    eta_plot = plottools.pcolorcells(fgout1.X,fgout1.Y,eta,
                                     cmap=geoplot.tsunami_colormap)
    clim(-4,4)
    cb = colorbar(eta_plot, extend='both', shrink=0.5)
    cb.set_label('meters')
    title_text = title('Surface at time %s\n%s' \
            % (timedelta(seconds=fgout1.t), event))

    ax.set_aspect(1./cos(ylat*pi/180.))
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    ax.set_xlim(plot_extent[:2])
    ax.set_ylim(plot_extent[2:])

            
            
    def update(fgframeno):
        """
        Update an exisiting plot with solution from fgout frame fgframeno.
        """
        
        fgout = fgout_grid.read_frame(fgframeno)
        print('Updating plot at time %s' % timedelta(seconds=fgout.t))
            
        # reset title to current time:
        title_text.set_text('Surface at time %s\n%s' \
                    % (timedelta(seconds=fgout.t),event))

        # reset surface eta to current state:
        eta = ma.masked_where(fgout.h<0.001, fgout.eta)
        eta_plot.set_array(eta.T.flatten())


    def plot_fgframe(fgframeno):
        """
        Convenience function for plotting one frame.
        But if you use this function in IPython and then try to make the animation,
        it may get into an infinite loop (not sure why).  Close the figure to abort.
        """
        update(fgframeno)
                    

    def make_anim2():
        print('Making anim...')
        anim = animation.FuncAnimation(fig, update,
                                       frames=fgframes, 
                                       interval=200, blit=False)
        return anim


    #anim = make_anim2()
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   interval=200, blit=False)
    
    # Output files:
    name = '%s_%s_animation' % (location,event)

    fname_mp4 = os.path.join(plotdir, name + '.mp4')
    fname_html = None
    #fname_html = os.path.join(plotdir, name + '.html')
    
    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        fname_html = name + '.html'
        animation_tools.make_html(anim, file_name=fname_html, title=name)
        
        
if __name__ == '__main__':
    
    event = os.path.split(os.getcwd())[-1]

    outdir = os.path.abspath('./_output')
    plotdir = os.path.abspath('./_plots')

    make_anim(outdir,plotdir,location,event)
    
