"""
Make an mp4 animation of fgout grid results. 
Comparison of time-dependent and instantaneous ruptures

JUST STARTED - want just one colorbar

This is done in a way that makes the animation quickly and with minimum 
storage required, by making one plot and then defining an update function
that only changes the parts of the plot that change in each frame.
The tuple update_artists contains the list of Artists that must be changed
in update.  Modify this as needed.

"""

from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from datetime import timedelta

from clawpack.geoclaw import fgout_tools
    
event = os.path.split(os.getcwd())[-1]
event_instant = event + '_instant'

fgno = 1  # which fgout grid

outdir = '_output'
outdir_instant = '../%s_instant/_output' % event
format = 'binary'  # format of fgout grid output

#fgframes = range(1,241)  # frames of fgout solution to use in animation
fgframes = range(1,121)  # 1 hour
#fgframes = range(1,61,15)

figsize = (10,8)

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 
fgout_grid_instant = fgout_tools.FGoutGrid(fgno, outdir_instant, format) 

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

fgout = fgout_grid.read_frame(fgframes[0])
fgout_instant = fgout_grid_instant.read_frame(fgframes[0])

#plot_extent = fgout.extent_edges
plot_extent = [-130,-122.5,38.5,50.5]

ylat = fgout.Y.mean()  # for aspect ratio of plots

fig = figure(figsize=figsize)
ax0 = axes([.1,.1,.35,.8])
ax1 = axes([.5,.1,.35,.8])
cax = axes([.88,.2,.02,.6])  # colorbar

axs = [ax0,ax1]
eta_plot = [0,0]

for k,ax in enumerate(axs):
    ax.set_xlim(plot_extent[:2])
    ax.set_ylim(plot_extent[2:])

    topo_plot = plottools.pcolorcells(fgout.X,fgout.Y,fgout.B, ax=ax,
                                      cmap=geoplot.land1_colormap)
    topo_plot.set_clim(0,2000)

    if k==0:
        eta = ma.masked_where(fgout.h<0.001, fgout.eta)
        eta_plot0 = plottools.pcolorcells(fgout.X,fgout.Y,eta,
                                     ax=ax,cmap=geoplot.tsunami_colormap)
        eta_plot0.set_clim(-4,4)
    else:
        eta = ma.masked_where(fgout_instant.h<0.001, fgout_instant.eta)
        eta_plot1 = plottools.pcolorcells(fgout_instant.X,fgout_instant.Y,eta,
                                     ax=ax,cmap=geoplot.tsunami_colormap)        
        eta_plot1.set_clim(-4,4)
        
    ax.set_aspect(1./cos(ylat*pi/180.))
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    ax.set_xlim(plot_extent[:2])
    ax.set_ylim(plot_extent[2:])

    

cb = colorbar(eta_plot1, cax=cax, extend='both')
cb.set_label('meters')

title_text0 = ax0.set_title('Surface at time %s\n%s' \
        % (timedelta(seconds=fgout.t), event))

title_text1 = ax1.set_title('Surface at time %s\n%s' \
        % (timedelta(seconds=fgout_instant.t), event_instant))


# The artists that will be updated for subsequent frames:
update_artists = (eta_plot0, title_text0, eta_plot1, title_text1)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
    
    fgout = fgout_grid.read_frame(fgframeno)
    fgout_instant = fgout_grid_instant.read_frame(fgframeno)
    print('Updating plot at time %s' % timedelta(seconds=fgout.t))
    
    # unpack update_artists (must agree with definition above):
    eta_plot0, title_text0, eta_plot1, title_text1 = update_artists
        
    # reset title to current time:
    title_text0.set_text('Surface at time %s\n%s' \
                % (timedelta(seconds=fgout.t),event))
    title_text1.set_text('Surface at time %s\n%s' \
                % (timedelta(seconds=fgout.t),event_instant))
                
    # reset surface eta to current state:
    eta = ma.masked_where(fgout.h<0.001, fgout.eta)
    eta_plot0.set_array(eta.T.flatten())
    
    eta = ma.masked_where(fgout_instant.h<0.001, fgout_instant.eta)
    eta_plot1.set_array(eta.T.flatten())
    
    update_artists = (eta_plot0, title_text0, eta_plot1, title_text1)
    return update_artists

def plot_fgframe(fgframeno):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':

    anim = make_anim()
    
    # Output files:
    event = os.path.split(os.getcwd())[-1]
    name = 'geoclaw_' + event + '_dual'

    fname_mp4 = name + '.mp4'
    fname_html = None # name + '.html'
    
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
    
    
    
    
