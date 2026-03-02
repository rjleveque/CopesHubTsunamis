"""
Plot fgout frames and transects
Specific to SBIT, for other sites it is necessary to change the axis
limits and transects, and provide a different background image.
"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

import matplotlib.pyplot as plt
import imageio_ffmpeg
plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()

from pylab import *
import os,sys,glob
from clawpack.visclaw import plottools, geoplot, gridtools
from clawpack.visclaw import animation_tools, colormaps
from matplotlib import animation, colors
from clawpack.geoclaw import fgout_tools
from datetime import timedelta
from clawpack.geoclaw import fgout_tools

def make_anim(outdir, plotdir, location, event):
    os.system('mkdir -p %s' % plotdir)
    print('Will take output from \n   %s and send plots to \n    %s' \
            % (outdir,plotdir))
    
    fgno = 1
    
    if 1:
        # all frames found in outdir:
        fgout_frames = glob.glob(os.path.join(outdir, \
                                              'fgout%s.t*' % str(fgno).zfill(4)))
        fgframes = range(1, len(fgout_frames)+1)
        nout = len(fgout_frames)
        print('Found %i fgout frames' % nout)
    
    if 0:
        nout = 301
        fgframes = range(1,nout+1)
        
    #fgframes = [1,nout]  # test on few frames
    
    print('Looking for output in ',outdir)
    
    output_format = 'binary32'
    
    # Instantiate object for reading fgout frames:
    fgout_grid1 = fgout_tools.FGoutGrid(fgno, outdir, output_format)
    fgout_grid1.read_fgout_grids_data()
    
    
    fgframes1 = array(fgframes)
    fgframes1 = [int(i) for i in fgframes1]  # convert from numpy.int64 to int
    
    print('fgframes1 = ',fgframes1)
    
    fgframe1 = fgframes1[0] # start with first frame
    
    fgout1 = fgout_grid1.read_frame(fgframe1)
    
    
    # ----------
    # Plotting:
    
    # Plot one frame of fgout data and define the Artists that will need to
    # be updated in subsequent frames:
    
    # Note that if you change what is plotted for each frame
    # (e.g. the transect locations) then you also need to change the
    # update function to update with the proper data each frame.
    
    #cmap_water = geoplot.tsunami_colormap  # opaque
    alpha = 0.8  # transparency
    cmap_water = colormaps.make_colormap({-0.6:[0,0,1,alpha],
                                           0.0:[0,1,1,alpha],
                                           0.6:[1,0,0,alpha]})
    
    fig = figure(1, figsize=(8,7))
    clf()
    
    
    plot_extent = fgout1.extent_edges
    
    ax = axes([.1,.4,.8,.5])

    if 0:
        ax.imshow(GE_image, extent=GE_extent)
        B = nan*fgout1.B
    else:
        B = fgout1.B
        
    #Here B started at time 0,so is B0
    B0 = copy(B)

    B_plot1 = ax.imshow(flipud(B.T), extent=fgout1.extent_edges,
           cmap=geoplot.land1_colormap)
           #cmap=geoplot.googleearth_transparent)
    
    B_plot1.set_clim(0,15)
    
    h_plus_B0 = fgout1.h + B0
    h_plus_B0_water = where(fgout1.h > 0.01,h_plus_B0,nan)
    h_water = where(fgout1.h > 0.01, fgout1.h, nan)
    eta_water = where(B0 > 0,h_water,h_plus_B0_water)
    eta_plot1 = ax.imshow(flipud(eta_water.T), extent=fgout1.extent_edges,
           cmap=cmap_water)
    
    climits = (-10,10)
    eta_plot1.set_clim(climits)
    axis(plot_extent)
    title_text = ax.set_title('%s\nSurface/Depth at time %s  (frame %i)' \
                % (event, timedelta(seconds=fgout1.t), fgframe1))
                
    ax.set_aspect(1/cos(46.73*pi/180))
    ticklabel_format(useOffset=False)
    
    
    if 1:
        cb = colorbar(eta_plot1, extend='both', shrink=0.7)
                      #orientation='vertical',anchor=(0,0))
                      #orientation='horizontal', anchor=(0.4,1))
        cb.set_label('meters')
        
    # Add transects to planview plot:
    # G3 is Gauge 3 at the intersection of the Casino SR105 and Tokeland Road
    yt1 = 46.72462963; Ttitle1 = 'G3'
    
    #west, east limits of the transect
    x1trans, x2trans = fgout1.extent_edges[:2]
    
    plot([x1trans,x2trans], [yt1,yt1],'k-',linewidth=0.8)
    text(-124.1,yt1+0.0005,'Transect %s' % Ttitle1)

    #plot([x1trans,x2trans], [yt2,yt2],'k-',linewidth=0.8)
    #text(x1trans-0.005,yt2+0.0005,'Transect 2 %s' % Ttitle2, fontsize=8)
    #plot([x1trans,x2trans], [yt3,yt3],'k-',linewidth=0.8)
    #text(x1trans-0.005,yt3+0.0005,'Transect 3 %s' % Ttitle3, fontsize=8)

    #end of plotting the plan view
    
    # =========
    # transects:
    # =========
    
    
    def extract_transect(fgout_soln,xtrans,ytrans):
    
        eta1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                       fgout_soln.eta, xtrans, ytrans)
        B1d = gridtools.grid_eval_2d(fgout_soln.X, fgout_soln.Y,
                                     fgout_soln.B, xtrans, ytrans)
        return B1d, eta1d
    
    def annotate_transect(axtrans):
        dxkm = 1
        dxlong = dxkm/(111.* cos(pi*46/180))
        axtrans.plot([-124.1,-124.1+dxlong], [-12,-12],'k')
        axtrans.text(-124.1+dxlong/2, -13, '%i km' % dxkm, \
                     ha='center',va='top')
    
        xG3 = -124.02037037; yG3= 46.72462963;
        axtrans.plot([xG3,xG3],[-20,0],'k--')
        axtrans.text(-124.015, -13, ' G3')
        axtrans.grid(True)
        axtrans.ticklabel_format(useOffset=False)
        axtrans.set_xlabel('longitude')
        axtrans.set_ylabel('meters')
        axtrans.set_xlim(x1trans,x2trans)
        xt = axtrans.get_xticks()
        #axtrans.set_xticks(xt,rotation=20)
        #axtrans.set_xticks(arange(x1trans,x2trans+1e-6,.01))
        
    ylimtr = (-20,20)  # ylimits for transect plots
    xtrans = linspace(x1trans, x2trans, 1000)  # x points on transects
    
    
    # Transect 1 (underneath plan view plot)
    
    y1trans, y2trans = 2*[yt1]; Ttitle = Ttitle1
    ytrans = linspace(y1trans, y2trans, 1000)
    
    #The four numbers in axes command are SW corner of axes, width, height in fractions of space
    #axtrans = axes([.55,.7,.35,.2])
    axtrans = axes([.1,.1,.8,.25])
    axtrans.set_title('Transect C at y = %.5f' % y1trans)
    
    axtrans.set_ylim(ylimtr)
    
    Btrans, etatrans = extract_transect(fgout1,xtrans,ytrans)
    #import pdb; pdb.set_trace()
    
    #Btrans, etatrans = Btrans1, etatrans1
    
    # filled regions:
    Bfill_plot = axtrans.fill_between(xtrans, Btrans-1e4, Btrans, 
                                      color=[.5,1,.5,1])
    etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, 
                                      color=[.5,.5,1,1])
    
    # surface and topo plots:
    etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
    Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')
    
    annotate_transect(axtrans)
    #end of  Transect 1

    figdummy,axdummy = subplots()
    
    def update(fgframeno):
        """
        Update an exisiting plot with solution from fgout frame fgframeno.
        """
            
        fgout1 = fgout_grid1.read_frame(fgframeno)
    
        print('Updating plot at time %s' % timedelta(seconds=fgout1.t))
            
        # reset title to current time:
        title_text.set_text('%s\nSurface/Depth at time %s  (frame %i)' \
                % (event,timedelta(seconds=fgout1.t), fgframeno))
    
        # reset eta and B in plan-view plots to current state:
    
        h_plus_B0 = fgout1.h + B0
        h_plus_B0_water = where(fgout1.h > 0.01, h_plus_B0, nan)
        h_water = where(fgout1.h > 0.01, fgout1.h, nan)
        eta_water = where(B0 > 0, h_water, h_plus_B0_water) # for zeta on GE

        eta_plot1.set_array(flipud(eta_water.T))

        if 0:
            eta_water = where(fgout1.h > 0.01, fgout1.eta, nan)
            h_water = where(fgout1.h > 0.01, fgout1.h, nan)
            eta_water = where(fgout1.B > 0, h_water, eta_water) # for zeta on GE
            eta_plot1.set_array(flipud(eta_water.T))
            #B_plot1.set_array(flipud(fgout1.B.T))
    
        # update transect:
        
        y1trans, y2trans = 2*[yt1]
    
        xtrans = linspace(x1trans, x2trans, 1000)
        ytrans = linspace(y1trans, y2trans, 1000)
        
        Btrans1, etatrans1 = extract_transect(fgout1,xtrans,ytrans)
        Btrans, etatrans = Btrans1, etatrans1
        Btrans_plot.set_data(xtrans,Btrans)
        etatrans_plot.set_data(xtrans,etatrans)
    
        #update the PolyCollections for fill_between plots:             
        dummy = axdummy.fill_between(xtrans, Btrans-1e4, Btrans, 
                                          color=[.5,1,.5,1])
        dp = dummy.get_paths()[0]
        dummy.remove()
        Bfill_plot.set_paths([dp.vertices])
    
        dummy = axdummy.fill_between(xtrans, Btrans, etatrans, 
                                          color=[.5,.5,1,1])
        dp = dummy.get_paths()[0]
        dummy.remove()
        etafill_plot.set_paths([dp.vertices])
    
    
    def plot_fgframe(fgframeno, save_png=False):
        """
        Convenience function for plotting one frame.
        But if you use this function in IPython and then try to make the animation,
        it may get into an infinite loop (not sure why).  Close the figure to abort.
        """
        update(fgframeno)
        
        if save_png:
            fname = 'fgout_frame%s.png' % str(fgframeno).zfill(4)
            savefig(fname, bbox_inches='tight')
            print('Created ',fname)
    
                    
    
    def make_anim():
        print('Making anim...')
        anim = animation.FuncAnimation(fig, update,
                                       frames=fgframes, 
                                       interval=200, blit=False)
        return anim
    
    
    anim = make_anim()

    # Output files:
    name = '%s_%s_animation' % (location,event)

    fname_mp4 = os.path.join(plotdir, name + '.mp4')
    #fname_html = os.path.join(plotdir, name + '.html')
    fname_html = None
    
    
    if fname_mp4:
        fps = 5
        print('Making mp4...')
        #writer = animation.writers['ffmpeg'](fps=fps)
        writer = animation.FFMpegWriter(fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        fname_html = name + '.html'
        animation_tools.make_html(anim, file_name=fname_html, title=name)

