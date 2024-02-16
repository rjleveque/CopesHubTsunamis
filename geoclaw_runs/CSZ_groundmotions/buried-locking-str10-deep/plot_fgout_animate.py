
from pylab import *
import os
from clawpack.visclaw import plottools, geoplot
from clawpack.visclaw import animation_tools
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta 

#fgno = 1
outdir = '_output'
output_format = 'binary'

# Instantiate object for reading fgout frames:
fgout_grid1 = fgout_tools.FGoutGrid(1, outdir, output_format)
fgout_grid2 = fgout_tools.FGoutGrid(2, outdir, output_format)

y0 = 46.90  # for transect
x1trans, x2trans = -124.5,-124.02


fgframes = range(1,122)
fgframe = fgframes[0] # start with first frame

fgout1 = fgout_grid1.read_frame(fgframe)
fgout2 = fgout_grid2.read_frame(fgframe)


#plot_extent = fgout.extent_edges
plot_extent = [-130,-122.5,38.5,50.5]

#fig,ax = subplots(num=1,, figsize=(13,8))
fig = figure(1, figsize=(13,8))
clf()
ax = axes([.1,.1,.4,.8])
B_plot1 = ax.imshow(flipud(fgout1.B.T), extent=fgout1.extent_edges,
       cmap=geoplot.land1_colormap)

B_plot1.set_clim(0,1500)

eta_water = where(fgout1.h > 0, fgout1.eta, nan)
eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
       cmap=geoplot.tsunami_colormap)

eta_plot1.set_clim(-3,3)
axis(plot_extent)
title_text = ax.set_title('Surface at time %s' % timedelta(seconds=fgout1.t))


if 1:
    cb = colorbar(eta_plot1, extend='both', shrink=0.4, 
                  orientation='vertical',anchor=(0,0))
                  #orientation='horizontal', anchor=(3,3))
    cb.set_label('meters')
    

# inset axes....
axins = ax.inset_axes([1.2,0.45,1.2,.7])
B_plot2 = axins.imshow(flipud(fgout2.B.T), extent=fgout2.extent_edges,
       cmap=geoplot.land_colors)
B_plot2.set_clim(0,100)
eta_water2 = where(fgout2.h > 0, fgout2.eta, nan)
eta_plot2 = axins.imshow(flipud(eta_water2.T), extent=fgout2.extent_edges,
       cmap=geoplot.tsunami_colormap)
eta_plot2.set_clim(-3,3)
axins.plot([x1trans,x2trans], [y0,y0], 'k', linewidth=0.7)
axins.text(-124.25, y0+0.005, 'Transect', fontsize=9)

# sub region of the original image
x1, x2, y1, y2 = -124.3,-123.75,46.8,47.1
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
#axins.set_xticklabels([])
axins.set_yticklabels([])
axins.secondary_yaxis('right')
axins.set_title('Grays Harbor')

ax.indicate_inset_zoom(axins, edgecolor="black")



# transect

axtrans = axes([.55,.1,.4,.3])
axtrans.set_title('Transect across Westport at latitude %.3f' % y0)

j1 = where(fgout1.y < y0)[0].max()
i1 = where(fgout1.x < -124.3)[0]
# where(fgout1.x < -124.3)[0]
#axtrans.plot(fgout1.x[i1], fgout1.B[i1,j1], 'g')
#eta_wet = where(fgout1.B[i1,j1]<0, fgout1.eta[i1,j1], nan)
eta_wet = fgout1.eta[i1,j1]
xtrans = fgout1.x[i1]
etatrans = eta_wet
Btrans = fgout1.B[i1,j1]
#axtrans.plot(fgout1.x[i1], eta_wet, 'b')
axtrans.set_xlim(x1trans,x2trans)
axtrans.set_ylim(-20,15)

j2 = where(fgout2.y < y0)[0].max()
i2 = where(fgout2.x >= -124.3)[0]
#eta_wet = where(fgout2.B[i2,j2]<0, fgout2.eta[i2,j2], nan)
eta_wet = fgout2.eta[i2,j2]

xtrans = hstack((xtrans, fgout2.x[i2]))
etatrans = hstack((etatrans,eta_wet))
Btrans = hstack((Btrans,fgout2.B[i2,j2]))

Bfill_plot = axtrans.fill_between(xtrans, Btrans-20, Btrans, color=[.5,1,.5,0])
etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, color=[.5,.5,1,0])


etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')


axtrans.grid(True)



# The artists that will be updated for subsequent frames:
update_artists = (B_plot1, eta_plot1, B_plot2, eta_plot2,
                  Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot,
                  title_text)
        
        
def update(fgframeno, *update_artists):
    """
    Update an exisiting plot with solution from fgout frame fgframeno.
    The artists in update_artists must have new data assigned.
    """
        
    fgout1 = fgout_grid1.read_frame(fgframeno)
    fgout2 = fgout_grid2.read_frame(fgframeno)

    print('Updating plot at time %s' % timedelta(seconds=fgout1.t))
    
    # unpack update_artists (must agree with definition above):
    B_plot1, eta_plot1, B_plot2, eta_plot2, \
          Bfill_plot, etafill_plot, Btrans_plot, etatrans_plot, \
          title_text = update_artists
        
    # reset title to current time:
    title_text.set_text('Surface at time %s' % timedelta(seconds=fgout1.t))

    # reset eta and B in plots to current state:
    
    eta_water = where(fgout1.h > 0, fgout1.eta, nan)
    eta_plot1.set_array(flipud(eta_water.T))
    B_plot1.set_array(flipud(fgout1.B.T))
            
    eta_water = where(fgout2.h > 0, fgout2.eta, nan)
    eta_plot2.set_array(flipud(eta_water.T))
    B_plot2.set_array(flipud(fgout2.B.T))

    #j1 = where(fgout1.y < y0)[0].max()
    #i1 = where(fgout1.x < -124.3)[0]
    #eta_wet = where(fgout1.B[i1,j1]<0, fgout1.eta[i1,j1], nan)
    eta_wet = fgout1.eta[i1,j1]
    
    xtrans = fgout1.x[i1]
    etatrans = eta_wet
    Btrans = fgout1.B[i1,j1]

    #j2 = where(fgout2.y < y0)[0].max()
    #i2 = where(fgout2.x >= -124.3)[0]
    #eta_wet = where(fgout2.B[i2,j2]<0, fgout2.eta[i2,j2], nan)
    eta_wet = fgout2.eta[i2,j2]
    
    xtrans = hstack((xtrans, fgout2.x[i2]))
    etatrans = hstack((etatrans,eta_wet))
    Btrans = hstack((Btrans,fgout2.B[i2,j2]))
    
    #Bfill_plot.set_data(xtrans, Btrans-20, Btrans)
    #etafill_plot.set_data(xtrans, B, etatrans)
    #Bfill_plot = axtrans.fill_between(xtrans, Btrans-20, Btrans, color=[.5,1,.5])
    #etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, color=[.5,.5,1])

    Btrans_plot.set_data(xtrans,Btrans)
    etatrans_plot.set_data(xtrans,etatrans)
    
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
                                   frames=fgframes, 
                                   fargs=update_artists,
                                   interval=200, blit=True)
    return anim

if __name__ == '__main__':
    
    anim = make_anim()
    
    # Output files:
    name = 'fgout_animation'

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