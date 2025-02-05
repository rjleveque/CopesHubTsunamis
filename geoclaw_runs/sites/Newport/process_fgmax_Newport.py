
"""
Process fgmax grid results and make plots of:
    preseismic topography B0 and postseismic B
    maximum surface elevation offshore
    maximum flow depth onshore (based on where B0 > 0)
    maximum flow speed
#Also create a kmz file with plots to be viewed on Google Earth.

Before running this code, make sure the file specified by fname_B0 exists.
This file provides the pre-seismic topography at each fgmax point.
See README for instructions on creating it by doing a GeoClaw run with no dtopo.

Now includes code to rewrite some of the fgmax data as a netCDF file that is
smaller and easier for others to use.
"""


import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *

import os,sys,glob,zipfile,shutil
from clawpack.clawutil.data import ClawData
from clawpack.geoclaw import topotools, dtopotools
from clawpack.visclaw import colormaps, gridtools
import matplotlib as mpl
from matplotlib import colors
from clawpack.amrclaw import region_tools
from clawpack.visclaw import plottools
from clawpack.geoclaw import kmltools, fgmax_tools
import sys


save_figs = True             # make png files for figures?
close_figs = True            # close big figures after saving?

## Also see name==main where it is read in from params.py
#location = 'Newport'

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Set CHT enviornment variable to repository top")
    
use_force_dry = False
if use_force_dry:
    fname_force_dry = os.path.join(input_dir, 'force_dry_init.data')
    print('Using force_dry_init from ', fname_force_dry)


def load_fgmax(outdir,fgno,fname_B0):
    # Read fgmax data:
    fg = fgmax_tools.FGmaxGrid()
    fgmax_input_file_name = outdir + '/fgmax_grids.data'
    print('fgmax input file: \n  %s' % fgmax_input_file_name)
    fg.read_fgmax_grids_data(fgno=fgno, data_file=fgmax_input_file_name)

    # determine time interval used for fgmax:
    clawdata = ClawData()
    clawdata.read(os.path.join(outdir, 'claw.data'), force=True)
    try:
        t_hours = min(clawdata.tfinal, fg.tend_max) / 3600.
    except:
        t_hours = nan

    fg.read_output(outdir=outdir, indexing='xy')  # so array layout same as topofile

    #### Read pre-seismic B0 from special run with no dtopo specified

    if 1:
        topoB0 = topotools.Topography()
        topoB0.read(fname_B0, topo_type=3)
        B0 = topoB0.Z
        B0_masked = ma.masked_array(B0, fg.B.mask)
        fg.B0 = B0_masked
    else:
        fg.B0 = fg.B
        print('No subsidence or uplift in fgmax region')

    dB = fg.B - fg.B0
    print('Minimum/maximum dB in fgmax region: %.2f m, %.2f m' \
            % (dB.min(), dB.max()))
            
    return fg, t_hours

# colormaps used below:

# colormap for depth:
bounds_depth = array([1e-6,1,2,4,6,10,12])

cmap_depth = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],
                 [1,.7,.7], [1,.4,.4], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_depth.set_over(color=[1,0,1])

# Set color for land points without inundation to light green:
#cmap_depth.set_under(color=[.7,1,.7])
# Set color for land points without inundation to transparent if on image:
cmap_depth.set_under(color=[.7,1,.7,0])

norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

# colormap for eta:
bounds_eta = array([0,1,2,4,6,10,12])

cmap_eta = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],
                 [1,.7,.7], [1,.4,.4], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_eta.set_over(color=[1,0,1])

norm_eta = colors.BoundaryNorm(bounds_eta, cmap_eta.N)

# colormap for speed:
bounds_speed = np.array([1e-6,1,2,4,6,8,10,12])
cmap_speed = mpl.colors.ListedColormap([[.9,.9,1],[.6,.6,1],                
                [.3,.3,1],[0,0,1], [1,.8,.8],
                [1,.6,.6], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_speed.set_over(color=[1,0,1])

# Set color for land points without inundation to light green:
#cmap_speed.set_under(color=[.7,1,.7])
# Set color for land points without inundation to transparent if on image:
cmap_speed.set_under(color=[.7,1,.7,0])

norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)
    


def make_fgmax_plots(fg, fgmax_plotdir, run_name, t_hours, GE_image, GE_extent):

    #Here, GE_extent is the fgmax_extent
    fgmax_extent = GE_extent
    
    def savefigp(fname):
        global save_figs
        if save_figs:
            fullname = '%s/%s_%s' % (fgmax_plotdir, run_name, fname)
            savefig(fullname, bbox_inches='tight')
            print('Created ', fullname)
        else:
            print('save_figs = False')

    ylat = fg.y.mean()  # for aspect ratio of plots

    ##### Plot topography

    zmin = -10.
    zmax = 10.
    land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                         0.25:[0.0,1.0,0.0],
                                          0.5:[0.8,1.0,0.5],
                                          1.0:[0.8,0.5,0.2]})

    sea_cmap = colormaps.make_colormap({ 0.0:[0,0,1], 1.:[.8,.8,1]})

    cmap, norm = colormaps.add_colormaps((land_cmap, sea_cmap),
                                         data_limits=(zmin,zmax),
                                         data_break=0.)                                   

    def plotZ(Z, show_cb=True):
        pc = plottools.pcolorcells(fg.X, fg.Y, Z, cmap=cmap, norm=norm)  
        if show_cb:
            cb = colorbar(pc,shrink=0.5,extend='both')
            cb.set_label('meters')
        gca().set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20);
        
        
    figure(figsize=(12,8))
    subplot(121)
    plotZ(fg.B0, show_cb=True)
    title('GeoClaw B0 before quake')

    subplot(122)
    plotZ(fg.B, show_cb=True)
    dB = fg.B - fg.B0
    #title('GeoClaw B after quake\nAverage subsidence = %.2f m' % dB.mean())
    title('GeoClaw B after quake\ndB.min = %.2f, dB.max = %.2f' % \
            (dB.min(),dB.max()))
    tight_layout()
    savefigp('geoclaw_topo.png')

    if 1: #For Newport
        othercondition = None    #onshore determined only by fg.B0 > 0

    if 0: #For Seaside
          # For Seaside, consider harbor/rivers to be onshore when plotting zeta:
        othercondition = fg.X > -123.93

    onshore = logical_or(fg.B0 >  0., othercondition)

    if fg.force_dry_init is not None:
        onshore = logical_or(onshore, fg.force_dry_init)
    offshore = logical_not(onshore)

    fg.h_onshore = ma.masked_where(offshore, fg.h)

    # zeta = h where B0>0 or h+B0 where B0<0 shows depth or apparent change in
    # water level in river or harbor when viewed from shore:

    fg.zeta_onshore = where(logical_and(fg.B0 <= 0., othercondition),
                            fg.B0+fg.h, fg.h_onshore)
    fg.zeta_onshore = ma.masked_where(offshore, fg.zeta_onshore)

    # use B0 for continuity at shore:                                    
    fg.eta_offshore = ma.masked_where(onshore, fg.B0 + fg.h)


    # Plot maximum flow depth

    maxh_onshore = nanmax(fg.zeta_onshore)
    maxh_onshore_ft = maxh_onshore/0.3048

    figure(figsize=(8,8))
    imshow(GE_image, extent=GE_extent)
    pc = plottools.pcolorcells(fg.X, fg.Y, fg.zeta_onshore, cmap=cmap_depth, norm=norm_depth)
    cb = colorbar(pc, extend='max', shrink=0.7)
    cb.set_label('meters')
    #contour(fg.X, fg.Y, fg.B0, [0], colors='g')

    gca().set_aspect(1./cos(ylat*pi/180.))
    axis(fgmax_extent)
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    title('Maximum onshore flow depth h over %.2f hours\n' % t_hours \
            +'(h+B0 in harbor/rivers), max = %.2f meters ' % maxh_onshore)
    
    #For Newport
    # Add transects to planview plot:
    one_third = 1.0/(3.0*3600.)
    yt1 = 44.635; Ttitle1 = ' '
    yt2 = 44.615; Ttitle2 = '(Yaq. Bay)'
    yt3 = 44.6025; Ttitle3 = ' '
    x1trans, x2trans = GE_extent[0]+one_third, GE_extent[1]-one_third 
    plot([x1trans,x2trans], [yt1,yt1],'-', color='yellow', linewidth=1.2)
    text(x1trans,yt1+0.0005,'Transect 1 %s' % Ttitle1, fontsize=12,
         ha='left', color='yellow')
    plot([x1trans,x2trans], [yt2,yt2],'-', color='yellow', linewidth=1.2)
    text(x1trans,yt2+0.0005,'Transect 2 %s' % Ttitle2, fontsize=12,
         ha='left', color='yellow')
    plot([x1trans,x2trans], [yt3,yt3],'-', color='yellow', linewidth=1.2)
    text(x1trans,yt3+0.0005,'Transect 3 %s' % Ttitle3, fontsize=12,
         ha='left', color='yellow')
    savefigp('h_onshore.png')
    #savefigp('zeta_onshore.png')


    # plot max speed
    figure(figsize=(8,8))
    imshow(GE_image, extent=GE_extent)
    pc = plottools.pcolorcells(fg.X, fg.Y, fg.s, cmap=cmap_speed, norm=norm_speed)
    cb = colorbar(pc, extend='max', shrink=0.7)
    cb.set_label('m/s')
    #contour(fg.X, fg.Y, fg.B0, [0], colors='g')
    gca().set_aspect(1./cos(ylat*pi/180.))
    axis(fgmax_extent)
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    maxs = fg.s.max()
    title('Maximum speed over %.2f hours\n' % t_hours \
        + 'max = %.2f m/s ' % maxs)
    savefigp('speed.png')


    # Plot eta offshore if desired:

    if 1:
        figure(figsize=(8,8))
        imshow(GE_image, extent=GE_extent)
        pc = plottools.pcolorcells(fg.X, fg.Y, fg.eta_offshore, cmap=cmap_eta, norm=norm_eta)
        cb = colorbar(pc, extend='max', shrink=0.7)
        cb.set_label('meters')
        #contour(fg.X, fg.Y, fg.B0, [0], colors='g')
        gca().set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        title('Maximum offshore surface eta over %.2f hours' % t_hours)
        savefigp('eta_offshore.png')

    #import pdb; pdb.set_trace()

    #Use same ones from before to see how plots look
    x1trans, x2trans = GE_extent[0]+one_third, GE_extent[1]-one_third 

    def extract_transect(fgmax_soln,xtrans,ytrans):
        h1d = gridtools.grid_eval_2d(fgmax_soln.X.T, fgmax_soln.Y.T,
                                       fgmax_soln.h.T, xtrans, ytrans)
        B1d = gridtools.grid_eval_2d(fgmax_soln.X.T, fgmax_soln.Y.T,
                                     fgmax_soln.B0.T, xtrans, ytrans)
        eta1d = h1d + B1d
        return B1d, eta1d

    fig = figure(1, figsize=(7,8))
    clf()

    def annotate_transect(axtrans):
        dxkm = 1
        dxlong = dxkm/(111.* cos(pi*ylat/180))
        axtrans.plot([x1trans,x1trans+dxlong], [-12,-12],'k')
        axtrans.text(x1trans+dxlong/2, -13, '%i km' % dxkm, \
                     ha='center',va='top')

        axtrans.grid(True)
        axtrans.ticklabel_format(useOffset=False)
        axtrans.set_xlabel('longitude')
        axtrans.set_ylabel('meters')
        axtrans.set_xlim(x1trans,x2trans)
        xt = axtrans.get_xticks()
        #axtrans.set_xticks(xt,rotation=20)
        axtrans.set_xticks(arange(x1trans,x2trans+1e-6,.01))
        
    ylimtr = (-20,25)  # ylimits for transect plots
    xtrans = linspace(x1trans, x2trans, 1000)  # x points on transects


    # Transect 1 (top)

    y1trans, y2trans = 2*[yt1]; Ttitle = Ttitle1
    ytrans = linspace(y1trans, y2trans, 1000)
    Btrans, etatrans = extract_transect(fg,xtrans,ytrans)

    axtrans = axes([.1,.7,.8,.2])
    axtrans.set_title('Transect 1 at y = %.5f %s - max depth on original topo' \
                        % (y1trans,Ttitle1))
    axtrans.set_ylim(ylimtr)

    # filled regions:
    Bfill_plot = axtrans.fill_between(xtrans, Btrans-1e4, Btrans, 
                                      color=[.5,1,.5,1])
    etafill_plot = axtrans.fill_between(xtrans, Btrans, etatrans, 
                                      color=[.5,.5,1,1])

    # surface and topo plots:
    etatrans_plot, = axtrans.plot(xtrans, etatrans, 'b')
    Btrans_plot, = axtrans.plot(xtrans, Btrans, 'g')
    annotate_transect(axtrans)


    # Transect 2 (middle)

    y1trans, y2trans = 2*[yt2]; Ttitle = Ttitle2
    ytrans = linspace(y1trans, y2trans, 1000)
    axtrans2 = axes([.1,.4,.8,.2],sharex=axtrans)
    axtrans2.set_title('Transect 2 at y = %.5f %s - max depth on original topo' \
                        % (y1trans,Ttitle2))

    axtrans2.set_ylim(ylimtr)

    Btrans2, etatrans2 = extract_transect(fg,xtrans,ytrans)

    # filled regions:
    Bfill_plot2 = axtrans2.fill_between(xtrans, Btrans2-1e4, Btrans2, 
                                      color=[.5,1,.5,1])
    etafill_plot2 = axtrans2.fill_between(xtrans, Btrans2, etatrans2, 
                                      color=[.5,.5,1,1])
    # surface and topo plots:
    etatrans_plot2, = axtrans2.plot(xtrans, etatrans2, 'b')
    Btrans_plot2, = axtrans2.plot(xtrans, Btrans2, 'g')

    annotate_transect(axtrans2)

    # Transect 3 (bottom)

    y1trans, y2trans = 2*[yt3]; Ttitle = Ttitle3

    ytrans = linspace(y1trans, y2trans, 1000)

    axtrans3 = axes([.1,.1,.8,.2],sharex=axtrans)
    axtrans3.set_title('Transect 3 at y = %.5f %s - max depth on original topo' \
                        % (y1trans,Ttitle3))

    #axtrans.set_xlim(x1trans,x2trans)
    #axtrans.sharex(ax)

    axtrans3.set_ylim(ylimtr)

    Btrans3, etatrans3 = extract_transect(fg,xtrans,ytrans)

    # filled regions:
    Bfill_plot3 = axtrans3.fill_between(xtrans, Btrans3-1e4, Btrans3, 
                                      color=[.5,1,.5,1])
    etafill_plot3 = axtrans3.fill_between(xtrans, Btrans3, etatrans3, 
                                      color=[.5,.5,1,1])
                                      
    # surface and topo plots:
    etatrans_plot3, = axtrans3.plot(xtrans, etatrans3, 'b')
    Btrans_plot3, = axtrans3.plot(xtrans, Btrans3, 'g')

    annotate_transect(axtrans3)

    savefigp('transects.png')

def make_kmz_plots(fg, fgmax_plotdir, run_name):
    # ## Plots for Google Earth overlays
    # 
    # The new version of `kmltools` includes some tools to make png files
    # that display properly on Google Earth.  The png files have no axes
    # and have the dimension and dpi set properly so that there is an integer
    # number of pixels in each grid cell so cell edges are sharp when zooming in.
    # 
    # 
    # We make three png files and then make a kml file that can be used to open all three.

    if 1:

        kml_dir = fgmax_plotdir + '/kmlfiles'
        print('Will send kml file and plots to kml_dir = \n  ', kml_dir)
        os.system('mkdir -p %s' % kml_dir);

        #fg.x = fg.X[:,0]
        #fg.y = fg.Y[0,:]

        h_wet_onshore = ma.masked_where(fg.h_onshore==0., fg.zeta_onshore)
        print('fg.x, fg.y shapes: ',fg.x.shape, fg.y.shape)
        print('+++ h_wet_onshore.shape = ',h_wet_onshore.shape)
        png_filename=kml_dir+'/h_onshore_max_for_kml.png'
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                         h_wet_onshore,
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_depth, norm=norm_depth)
        if close_figs: close('all')


        png_filename=kml_dir+'/eta_offshore_max_for_kml.png'
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                         fg.eta_offshore,
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_eta, norm=norm_eta)
        if close_figs: close('all')



        speed = ma.masked_where(fg.h==0., fg.s)
        png_filename = '%s/speed_max_for_kml.png' % kml_dir
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y, speed, 
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_speed, norm=norm_speed)
        if close_figs: close('all')


        stays_dry = ma.masked_where(fg.h>0., fg.h)
        png_filename = '%s/stays_dry_for_kml.png' % kml_dir
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                         stays_dry, 
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_speed, norm=norm_speed)
        if close_figs: close('all')


        # ### Make colorbars for kml files


        kmltools.kml_build_colorbar('%s/colorbar_depth.png' % kml_dir, cmap_depth, 
                                   norm=norm_depth, label='meters', title='depth', extend='max')
        kmltools.kml_build_colorbar('%s/colorbar_speed.png' % kml_dir, cmap_speed, 
                                   norm=norm_speed, label='meters / second', title='speed', extend='max')
        kmltools.kml_build_colorbar('%s/colorbar_eta.png' % kml_dir, cmap_eta, 
                                   norm=norm_eta, label='meters', title='eta', extend='max')
        if close_figs: close('all')


        # ### Make the kml file to display these three png files
        # 
        # Then you can open `fgmax_results_kmlfiles/fgmax_results.kml` in Google Earth to view them.



        if 1:
            # include eta offshore:
            png_files=['h_onshore_max_for_kml.png', 'speed_max_for_kml.png','stays_dry_for_kml.png',
                       'eta_offshore_max_for_kml.png']
            png_names=['max depth (zeta) onshore','max speed','stays dry',
                       'eta_offshore']
            cb_files = ['colorbar_depth.png', 'colorbar_speed.png',
                        'colorbar_eta.png']
            cb_names = ['colorbar_depth', 'colorbar_speed',
                        'colorbar_eta']

        else:
            # without eta offshore:
            png_files=['h_onshore_max_for_kml.png', 'speed_max_for_kml.png','stays_dry_for_kml.png']
            png_names=['max depth (zeta) onshore','max speed','stays dry']
            cb_files = ['colorbar_depth.png', 'colorbar_speed.png']
            cb_names = ['colorbar_depth', 'colorbar_speed']
                    
        name = 'fgmax_%s' % run_name
        fname = os.path.join(kml_dir, name+'.kml')
        kmltools.png2kml(png_extent, png_files=png_files, png_names=png_names, 
                         name=name, fname=fname,
                         radio_style=False,
                         cb_files=cb_files, cb_names=cb_names)


        # ## Create .kmz file including all plots

        savedir = os.getcwd()
        os.chdir(kml_dir)
        files = glob.glob('*.kml') + glob.glob('*.png')
        print('kmz file will include:')
        for file in files:
            print('    %s' % os.path.split(file)[-1])

        fname_kmz = '%s_fgmax.kmz' % run_name
        with zipfile.ZipFile(fname_kmz, 'w') as zip:
            for file in files:
                zip.write(file) 
            
        path_kmz = os.path.join(fgmax_plotdir, fname_kmz)
        shutil.move(fname_kmz, path_kmz)
        print('Created %s' % os.path.abspath(path_kmz))
        shutil.rmtree(kml_dir)
        os.chdir(savedir)

def make_nc_input(fname_nc, fg, force=False, verbose=True):

    import netCDF4
    import time
        
    if os.path.isfile(fname_nc):
        if force and verbose:
            print('Overwriting ', fname_nc)
        elif not force:
            print('*** netCDF file already exists, \n'\
                + '*** NOT overwriting '\
                + '--- use force==True to overwrite' )
            return -1
    
    with netCDF4.Dataset(fname_nc, 'w') as rootgrp:

        rootgrp.description = "fgmax data for " + fg.id
        rootgrp.history = "Created with input data " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
            
        if fg.X is not None:
            x = fg.X[0,:]
            lon = rootgrp.createDimension('lon', len(x))
            longitudes = rootgrp.createVariable('lon','f8',('lon',))
            longitudes[:] = x
            longitudes.units = 'degrees_east'
        else:
            if verbose: print('fg.X is None, not adding x')
            
        if fg.Y is not None:
            y = fg.Y[:,0]
            lat = rootgrp.createDimension('lat', len(y))
            latitudes = rootgrp.createVariable('lat','f8',('lat',))
            latitudes[:] = y
            latitudes.units = 'degrees_north'
        else:
            if verbose: print('fg.Y is None, not adding y')
            
        if fg.fgmax_point is not None:
            fgmax_point_var = \
                rootgrp.createVariable('fgmax_point','u1',('lat','lon',))
            fgmax_point_var[:,:] = fg.fgmax_point
        else:
            if verbose: print('fg.fgmax_point is None, not adding')
            
        if fg.force_dry_init is not None:
            force_dry_init = \
                rootgrp.createVariable('force_dry_init','u1',('lat','lon',))
            force_dry_init[:,:] = fg.force_dry_init
        else:
            if verbose: print('fg.force_dry_init is None, not adding')  

        print('Created %s' % fname_nc)            
        if verbose:
            print('History:  ', rootgrp.history) 
        return 0     
        
def write_nc_output(fname_nc, fg, new=False, force=False, 
                    outdir='Unknown', verbose=True):

    from clawpack.clawutil.data import ClawData 
    import netCDF4
    import time
    
    fv = -9999.   # fill_value for netcdf4
    
    if new:
        # first create a new .nc file with X,Y,fgmax_point,force_dry_init:
        result = make_nc_input(fname_nc, fg, force=force, verbose=verbose)
        if result == -1:
            print('*** make_nc_input failed, not appending output')
            return        
        
    if outdir == 'Unknown':
        # Cannot determine tfinal or run_finished time
        tfinal = fv
        run_finished = 'Unknown'
    else:
        claw = ClawData()
        claw.read(outdir+'/claw.data', force=True)

        try:
            if claw.output_style==1:
                tfinal = claw.tfinal
            elif claw.output_style==2:
                tfinal = array(claw.output_times).max()
        except:
            tfinal = fv
        
        try:
            mtime = os.path.getmtime(outdir+'/timing.txt')
            run_finished = time.ctime(mtime) 
        except:
            run_finished = 'Unknown'
            
    # add fgmax output results to existing file
    print(os.getcwd())
    with netCDF4.Dataset(fname_nc, 'a') as rootgrp:
        if verbose:
            print('Appending data from fg to nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('        fg.id: ', fg.id)
        
        h = rootgrp.variables.get('h', None)
        if (h is not None) and (not force):
            print('*** netCDF file already contains output,\n'\
                + '*** NOT overwriting '\
                + '--- use force==True to overwrite' )
            return
                
        x = array(rootgrp.variables['lon'])
        y = array(rootgrp.variables['lat'])
        X,Y = meshgrid(x,y)
        try:
            fgmax_point = array(rootgrp.variables['fgmax_point'])
        except:
            fgmax_point = None
        bounding_box = [x.min(),x.max(),y.min(),y.max()]
        
        dx = x[1]-x[0]
        Xclose = allclose(fg.X, X, atol=0.1*dx)
        Yclose = allclose(fg.Y, Y, atol=0.1*dx)
        
        if (fg.X.shape != X.shape):
            # for now raise an exception, might want to extent to allow
            # filling only part of input arrays
            print('*** Mismatch of fg with data in nc file:')
            print('fg.X.shape = ',fg.X.shape)
            print('nc  X.shape = ',X.shape)
            print('fg.bounding_box = ',fg.bounding_box())
            print('nc  bounding_box = ',bounding_box)
            raise ValueError('*** Mismatch of fg with data in nc file')
    
        Xclose = allclose(fg.X, X, atol=0.1*dx)
        Yclose = allclose(fg.Y, Y, atol=0.1*dx)
        if (not (Xclose and Yclose)):
            raise ValueError('*** Mismatch of fg.X or fg.Y with data in nc file')
            

        rootgrp.history += "Added output " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
        
        rootgrp.tfinal = tfinal
        rootgrp.outdir = os.path.abspath(outdir)
        rootgrp.run_finished = run_finished
        
        #fgmax_point = rootgrp.variables.get('fgmax_point', None)

        if fg.dz is not None:
            try:
                dz = rootgrp.variables['dz']
            except:
                dz = rootgrp.createVariable('dz','f4',('lat','lon',),
                                            fill_value=fv)
            dz[:,:] = fg.dz
            dz.units = 'meters'
            if verbose: print('    Adding fg.dz to nc file')
        else:
            if verbose: print('fg.dz is None, not adding')

        if fg.B0 is not None:
            try:
                B0 = rootgrp.variables['B0']
            except:
                B0 = rootgrp.createVariable('B0','f4',('lat','lon',),
                                            fill_value=fv)
            B0[:,:] = fg.B0
            B0.units = 'meters'
            if verbose: print('    Adding fg.B0 to nc file')
        else:
            if verbose: print('fg.B0 is None, not adding')
            
        if fg.B is not None:
            try:
                B = rootgrp.variables['B']
            except:
                B = rootgrp.createVariable('B','f4',('lat','lon',),
                                            fill_value=fv)
            B[:,:] = fg.B
            B.units = 'meters'
            if verbose: print('    Adding fg.B to nc file')
        else:
            if verbose: print('fg.B is None, not adding')
                        
        if fg.h is not None:
            try:
                h = rootgrp.variables['h']
            except:
                h = rootgrp.createVariable('h','f4',('lat','lon',),
                                            fill_value=fv)
            h[:,:] = fg.h
            h.units = 'meters'
            if verbose: print('    Adding fg.h to nc file')
        else:
            if verbose: print('fg.h is None, not adding')
            
        if fg.s is not None:        
            try:
                s = rootgrp.variables['s']
            except:
                s = rootgrp.createVariable('s','f4',('lat','lon',),
                                            fill_value=fv)
            s[:,:] = fg.s
            s.units = 'meters/second'
            if verbose: print('    Adding fg.s to nc file')
        else:
            if verbose: print('fg.s is None, not adding')
            
        if fg.hss is not None:        
            try:
                hss = rootgrp.variables['hss']
            except:
                hss = rootgrp.createVariable('hss','f4',('lat','lon',),
                                            fill_value=fv)
            hss[:,:] = fg.hss
            hss.units = 'meters^3/sec^2'
            if verbose: print('    Adding fg.hss to nc file')
        else:
            if verbose: print('fg.hss is None, not adding')
            
        if fg.hmin is not None:        
            try:
                hmin = rootgrp.variables['hmin']
            except:
                hmin = rootgrp.createVariable('hmin','f4',('lat','lon',),
                                            fill_value=fv)
            # negate hmin so that it is minimum flow depth min(h):
            hmin[:,:] = -fg.hmin
            hmin.units = 'meters'
            if verbose: print('    Adding fg.hmin to nc file')
        else:
            if verbose: print('fg.hmin is None, not adding')
            
        if fg.arrival_time is not None:        
            try:
                arrival_time = rootgrp.variables['arrival_time']
            except:
                arrival_time = rootgrp.createVariable('arrival_time','f4',('lat','lon',),
                                            fill_value=fv)
            arrival_time[:,:] = fg.arrival_time
            arrival_time.units = 'seconds'
            if verbose: print('    Adding fg.arrival_time to nc file')
        else:
            if verbose: print('fg.arrival_time is None, not adding')
            
        print('Created %s' % fname_nc)
        if verbose:
            print('History:  ', rootgrp.history)
            print('\nMetadata:')
            print('  outdir:  ', rootgrp.outdir)
            print('  run_finished:  ', rootgrp.run_finished)
            print('  tfinal:  ', rootgrp.tfinal)

def read_nc(fname_nc, verbose=True):

    import netCDF4
    import time
    import os
    from numpy import ma
    from clawpack.geoclaw import fgmax_tools

    print('Reading fgmax data from %s' % fname_nc)
    
    def get_as_array(var, fgvar=None):
        if fgvar is None:
            fgvar = var
        a = rootgrp.variables.get(var, None)
        if a is not None:
            if verbose: print('    Loaded %s as fg.%s' % (var,fgvar))
            return array(a)
        else:
            if verbose: print('    Did not find %s for fg.%s' \
                                % (var,fgvar))
            return None
                    
    fg = fgmax_tools.FGmaxGrid()

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        if verbose:
            print('Reading data to fg from nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('History:  ', rootgrp.history)


                
        x = get_as_array('lon','x')
        y = get_as_array('lat','y')
        
        if (x is None) or (y is None):
            print('*** Could not create grid')
            return None
            
        X,Y = meshgrid(x,y)
        fg.X = X
        fg.Y = Y
        if verbose:
            print('    Constructed fg.X and fg.Y')
        
        # arrays defined everywhere
        fg.dz = get_as_array('dz')
        fg.force_dry_init = get_as_array('force_dry_init')
        fg.fgmax_point = get_as_array('fgmax_point') 
        
        if fg.fgmax_point is not None:
            # arrays defined only at fgmax points: return as masked arrays:
            fgmask = 1 - fg.fgmax_point  # mask points that are not fgmax pts
            fg.B = ma.masked_array(get_as_array('B'), fgmask)
            fg.B0 = ma.masked_array(get_as_array('B0'), fgmask)
            fg.h = ma.masked_array(get_as_array('h'), fgmask)
            fg.s = ma.masked_array(get_as_array('s'), fgmask)
            fg.hss = ma.masked_array(get_as_array('hss'), fgmask)
            fg.hmin = ma.masked_array(get_as_array('hmin'), fgmask)
            fg.arrival_time = ma.masked_array(get_as_array('arrival_time'), fgmask)
        else:
            fg.B = get_as_array('B')
            fg.B0 = get_as_array('B0')
            fg.h = get_as_array('h')
            fg.s = get_as_array('s')
            fg.hss = get_as_array('hss')
            fg.hmin = get_as_array('hmin')
            fg.arrival_time = get_as_array('arrival_time')
            
    if verbose:
        print('Returning FGmaxGrid object')
    return fg

if __name__== '__main__':
    
    sys.path.insert(0,'.')
    from params import event, location

    #import fgmax_tools  # uses local version with transposed arrays
                        # should appear in v5.10.0
                        

    #### Note: Are running from say buried-deep directory, where the _output is and where we want _plots
    outdir = os.path.abspath('./_output')
    #outdir = '/Users/rjl/scratch/CHT_runs/sites/seaside/multirun_tests/geoclaw_outputs/_output_buried-random-str10-shallow'
    plotdir = os.path.abspath('./_plots')
    os.system('mkdir -p %s' % plotdir)
    graphics_dir = os.path.join(CHT, 'geoclaw_runs/sites/Newport')
    B0_dir = os.path.join(CHT, 'geoclaw_runs/sites/Newport')

    ###########
    #Choose the fgnos to use
    fgnos = [1,2]
    image_names = [graphics_dir + '/Newport_fgmax0001GE.jpg',graphics_dir + '/Newport_fgmax0002GE.jpg']
    GE_extents = [[-124.1,-124.000092593,44.565,44.641944445], [-123.999907407,-123.92,44.565,44.641944445]]
    fnames_B0 = [B0_dir + '/fgmax0001_13s_B0.asc',B0_dir + '/fgmax0002_13s_B0.asc']
    ##########
        
    for fgno in fgnos:
        GE_image = imread(image_names[fgno-1])
        GE_extent = GE_extents[fgno-1]
        fname_B0 = fnames_B0[fgno-1]

        print('Will read fgmax results from outdir = \n  ', outdir)
        fgmax_plotdir = plotdir + '/fgmax' + str(fgno)
        print('Will send plots to fgmax_plotdir = \n  ', fgmax_plotdir)
        os.system('mkdir -p %s' % fgmax_plotdir);

        run_name = '%s_%s' % (location,event)
    
        fg, t_hours = load_fgmax(outdir,fgno,fname_B0)
        make_fgmax_plots(fg, fgmax_plotdir, run_name, t_hours, GE_image, GE_extent)

        if 1:
            make_kmz_plots(fg, fgmax_plotdir, run_name)

        close('all')

        if 0:
            fname_nc = '%s_fgmax' + 'str(fgno)' +'.nc' % run_name
            write_nc_output(fname_nc, fg, new=True, force=True, 
                        outdir=outdir, verbose=True)
                        
            fg2 = read_nc(fname_nc, verbose=True)  # test reading it back in
            print('max abs(B-B0) = %.2f' % abs(fg2.B-fg2.B0).max())
        
