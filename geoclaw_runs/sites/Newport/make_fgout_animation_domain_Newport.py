"""
Plot fgout frames  (adapted from Nu'u).
"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os, glob
from clawpack.visclaw import plottools, geoplot, gridtools
from clawpack.visclaw import animation_tools, colormaps
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 
from scipy.interpolate import interp1d

event = os.path.split(os.getcwd())[-1]
outdir = '_output'
# Output files:
name = 'animation_offshore_%s' % event

qmap = 'geoclaw'
if 'bouss' in event:
    qmap = 'geoclaw-bouss'


print('Looking for output in ',outdir)

output_format = 'binary32'


# Instantiate object for reading fgout frames:
fgout_grid1 = fgout_tools.FGoutGrid(1, outdir, output_format)
                                    #qmap=qmap)
fgout_grid1.read_fgout_grids_data()

fgout_grid2 = fgout_tools.FGoutGrid(2, outdir, output_format)
                                    #qmap=qmap)
fgout_grid2.read_fgout_grids_data()


# transect:

# Westport
#x1trans = -126.5; y1trans = 46.8
#x2trans = -124.0; y2trans = 46.8

# Seaside:
#x1trans = -126.; y1trans = 45.9931
#x2trans = -123.9; y2trans = 45.9931
 
# Newport: Use Yaquina Bay transect
one_sixth = 1.0/(6.0*3600)
west = -124.145 - one_sixth; east = -124.055 - one_sixth;
x1trans = west; y1trans = 44.615
x2trans = east; y2trans = 44.615
xtrans = linspace(x1trans, x2trans, 1000)
ytrans = linspace(y1trans, y2trans, 1000)


# fgout frames to include in animation:

if 0:
    nframes = fgout_grid1.nout
    err_msg = '*** expected the same fgout.nout for grids 1 and 2'
    if (fgout_grid2.nout < nframes):
        raise ValueError(err_msg)
    elif (fgout_grid2.nout > nframes):
        print(err_msg)
    
if 1:
    # all frames found in outdir:
    fgno = 1
    fgout_frames = glob.glob(os.path.join(outdir, \
                                          'fgout%s.t*' % str(fgno).zfill(4)))
    nframes = len(fgout_frames)
    print('Found %i fgout frames for fgno=%i' % (nframes,fgno))
        
fgframes1 = fgframes2 = range(1,nframes+1)
#fgframes1 = fgframes2 = range(1,104,1) # test

fgframe1 = fgframes1[0] # start with first frame
fgframe2 = fgframes2[0] # for fgout grid 2 

fgout1 = fgout_grid1.read_frame(fgframe1) # fgout grid offshore region
fgout2 = fgout_grid2.read_frame(fgframe2) # finer fgout grid for zoom


# ----------
# Plotting:

# Plot one frame of fgout data and define the Artists that will need to
# be updated in subsequent frames:

y0 = 46.  # latitude for aspect ratio

fig = figure(1, figsize=(10,8))
clf()

# =========
# Ocean:
# =========

#plot_extent = fgout.extent_edges
#plot_extent = [-130, -122, 38.5, 50.]  # full domain for offshore_gauges
#plot_extent = [-128, -122, 44.5, 47.5]  # for Seaside inundation
plot_extent = [-126, -122, 43.5, 46.5]  # for Newport inundation


ax = axes([.1,.3,.3,.45])
B_plot1 = ax.imshow(flipud(fgout1.B.T), extent=fgout1.extent_edges,
       cmap=geoplot.land1_colormap)

B_plot1.set_clim(0,1500)

eta_water = where(fgout1.h > 0, fgout1.eta, nan)
eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
       cmap=geoplot.tsunami_colormap)
ax.set_aspect(1/cos(y0*pi/180))

eta_plot1.set_clim(-5,5)
axis(plot_extent)
#title_text = ax.set_title('Event %s \ntime %.1f minutes' \
#        % (event,(fgout1.t/60.)))
title_text = ax.set_title('Surface at time %s\n%s' \
    % (timedelta(seconds=fgout1.t), event))


if 1:
    cb = colorbar(eta_plot1, extend='both', shrink=0.8, 
                  #orientation='vertical',anchor=(0,0))
                  orientation='horizontal', anchor=(0.5,0.8))
    cb.set_label('meters')
    
# =========
# Transect region
# =========


#onecolor = [0.3,0.3,1]
onecolor = [0.7,0.7,1]
cmap_onecolor = colormaps.make_colormap({0.:onecolor, 1.:onecolor})

# inset axes....
axins = ax.inset_axes([0.9,0.6,2.2,1.0])
#axins.imshow(GE_image, extent=GE_extent)
#axins.imshow(hillshade_image, extent=hillshade_extent)

B = ma.masked_where(fgout2.h < 0.01, fgout2.B)
#B = nan*fgout2.B  # don't show topo if using GE_image
B_plot2 = axins.imshow(flipud(B.T), extent=fgout2.extent_edges,
       cmap=geoplot.googleearth_transparent)
B_plot2.set_clim(0,8)

dry_depth = 0.02  # gdepth
eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
h_water2 = where(fgout2.h > dry_depth, fgout2.h, nan)
zeta_water2 = where(fgout2.B>0, h_water2, eta_water2)
eta_plot2 = axins.imshow(flipud(zeta_water2.T), extent=fgout2.extent_edges,
       cmap=geoplot.tsunami_colormap)
eta_plot2.set_clim(-5,5)
axins.plot([x1trans,x2trans], [y1trans,y2trans], 'k', linewidth=0.9)
axins.text(x1trans+0.0025, y2trans+0.005, 'Transect',
           ha='left',va='bottom',color='k', fontsize=10)
#xticks(rotation=20)

# sub region of the original image
#x1, x2, y1, y2 = [-126.5,-124.1,46.587,47.227]
#x1, x2, y1, y2 = [-126.5,-124.,46.5,47.5]  # Grays Harbor
#x1, x2, y1, y2 = [-126.0, -123.9, 45.7, 46.4]  # Seaside

west = fgout2.extent_edges[0]; east = fgout2.extent_edges[1];
south = fgout2.extent_edges[2]; north = fgout2.extent_edges[3];
x1, x2, y1, y2 = [west, east, south, north]  # Newport

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
#axins.set_xticks(arange(-156.182,-156.175,.002))
axins.ticklabel_format(useOffset=False)
axins.set_aspect(1/cos(y0*pi/180))
axins.tick_params(axis='x',rotation=20)

axins.set_yticklabels([])
axins_yaxis2 = axins.secondary_yaxis('right')
axins_yaxis2.ticklabel_format(useOffset=False)
axins.set_title("Offshore Seaside")
axins.contour(fgout2.X,fgout2.Y,fgout2.B,[1],colors='g')

# show where this region is on the Ocean plot:
ax.indicate_inset_zoom(axins, edgecolor="black")

    

# =========
# transect:
# =========


def extract_transect(fgout_soln):

    eta1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                   fgout_soln.eta, xtrans, ytrans)
    B1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                 fgout_soln.B, xtrans, ytrans)
    return B1d, eta1d

                                 

axtrans = axes([.5,.15,.45,.3])
axtrans.set_title('Transect at y = %.3f' % y2trans)

axtrans.set_xlim(x1trans,x2trans)
axtrans.set_ylim(-10,30)
axtrans.tick_params(axis='x',rotation=20)

Btrans1, etatrans1 = extract_transect(fgout1)
Btrans2, etatrans2 = extract_transect(fgout2)

Btrans, etatrans = Btrans2, etatrans2


# filled regions:
Bfill_plot = axtrans.fill_between(xtrans, Btrans-1e4, Btrans, 
                                  color=[.5,1,.5,1])
etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, 
                                  color=[.5,.5,1,1])

# surface and topo plots:
etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')


axtrans.grid(True)
axtrans.ticklabel_format(useOffset=False)



# The artists that will be updated for subsequent frames:
update_artists = (B_plot1, eta_plot1, B_plot2, eta_plot2,
                  Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,
                  title_text)
        
figdummy,axdummy = subplots()

def update(k, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
        
    #fgout1 = fgout_grid1.read_frame(fgframeno)
    #fgout2 = fgout_grid2.read_frame(2*fgframeno-1)

    fgout1 = fgout_grid1.read_frame(fgframes1[k])
    fgout2 = fgout_grid2.read_frame(fgframes2[k])

    print('Updating plot at time %s' % timedelta(seconds=fgout1.t))
    
    # unpack update_artists (must agree with definition above):
    B_plot1, eta_plot1, B_plot2, eta_plot2, \
          Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot, \
          title_text = update_artists
        
    # reset title to current time:
    title_text.set_text('Surface at time %s\n%s' \
        % (timedelta(seconds=fgout1.t), event))
    #title_text.set_text('Surface at time %s' % timedelta(seconds=fgout1.t))

    # reset eta and B in plan-view plots to current state:

    eta_water = where(fgout1.h > 0, fgout1.eta, nan)
    eta_plot1.set_array(flipud(eta_water.T))
    B_plot1.set_array(flipud(fgout1.B.T))

    eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
    #h_water2 = where(fgout2.h > 0.02, fgout2.h, nan)
    h_water2 = where(fgout2.h > dry_depth, fgout2.h, nan)
    zeta_water2 = where(fgout2.B>0, h_water2, eta_water2)
    #eta_water = where(fgout2.X > -156.184, eta_water, nan)
    eta_plot2.set_array(flipud(zeta_water2.T))
    #B_plot2.set_array(flipud(fgout2.B.T))
    

    # update transects:
    
    Btrans1, etatrans1 = extract_transect(fgout1)
    Btrans2, etatrans2 = extract_transect(fgout2)

    Btrans, etatrans = Btrans2, etatrans2

    Btrans_plot.set_data(xtrans,Btrans)
    etatrans_plot.set_data(xtrans,etatrans)
    

    #update the PolyCollections for fill_between plots:             
    dummy = axdummy.fill_between(xtrans, Btrans-1e4, Btrans, 
                                      color=[.5,1,.5,1])
    #import pdb; pdb.set_trace()
    try:
        dp = dummy.get_paths()[0]
    except:
        print('could not get_paths for transect plot')
        dp = None
        
    dummy.remove()
    if dp:
        Bfill_plot.set_paths([dp.vertices])

    dummy = axdummy.fill_between(xtrans, Btrans, etatrans, 
                                      color=[.5,.5,1,1])
    try:
        dp = dummy.get_paths()[0]
    except:
        print('could not get_paths for transect plot')
        dp = None
        
    dummy.remove()
    if dp:
        etafill_plot.set_paths([dp.vertices])
    
    #fig.canvas.draw_idle()
    
    update_artists = (B_plot1, eta_plot1, B_plot2, eta_plot2,
                      Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,
                      title_text)
    return update_artists

def plot_fgframe(fgframeno, save_png=False):
    """
    Convenience function for plotting one frame.
    But if you use this function in IPython and then try to make the animation,
    it may get into an infinite loop (not sure why).  Close the figure to abort.
    """
    update(fgframeno, *update_artists)
    
    if save_png:
        fname = 'fgout_frame%s.png' % str(fgframeno).zfill(4)
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

                

def make_anim():
    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                                   frames=range(len(fgframes1)), 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':
    
    anim = make_anim()
    
    rundir = os.getcwd()
    if 'mmfs1' in rundir:
        # on hyak:
        scrdir = rundir.replace('mmfs1/home','gscratch/tsunami')
        outdir = os.path.join(scrdir, '_output')
    plotdir = outdir.replace('_output','_plots')

    if 1:
        #fgout_plotdir = plotdir + '/animations'
        fgout_plotdir = plotdir
        os.system('mkdir -p %s' % fgout_plotdir)
    else:
        fgout_plotdir = '.'


    fname_mp4 = os.path.join(fgout_plotdir, name + '.mp4')
    #fname_html = os.path.join(fgout_plotdir, name + '.html')
    fname_html = None
    
    
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
