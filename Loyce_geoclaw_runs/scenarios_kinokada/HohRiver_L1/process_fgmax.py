
"""
Making process_fgmax.py into a function that can be used in setplot_case.py

### Using 47.74 for the latitude in the gca statements for Hoh
### same as in setplot_case.py -- check what using for time frames there

The function fgreport takes the following seven arguments as below
process_fgmax.fgreport(rundir,outdir, plotdir, location, event, dtopofile, run_name)

"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os,glob,zipfile,shutil
from clawpack.geoclaw import topotools, dtopotools
from clawpack.visclaw import colormaps
import matplotlib as mpl
from matplotlib import colors
from matplotlib import colors
from clawpack.amrclaw import region_tools
from clawpack.visclaw import plottools
from clawpack.geoclaw import kmltools
from clawpack.geoclaw import fgmax_tools

root_dir = os.environ['CHT']
rundir = os.getcwd()
print ('root_dir printed by process_fgmax was: ',root_dir)
print ('rundir printed by process_fgmax was: ',rundir)
print (' ')

def fgreport(fgno,rundir,outdir,plotdir,location,event,dtopofile,run_name):

    input_dir = os.path.join(rundir, 'input_files')
    fgmax_plotdir = plotdir + '/_other_figures'
    if not os.path.isdir(fgmax_plotdir):
        os.mkdir(fgmax_plotdir)

    save_figs = True             # make png files for figures?
    close_figs = True            # close big figures after saving?

    print ('-----')
    print('fgno is: ',fgno)
    print('Will read fgmax results from outdir = \n  ', outdir)
    print('Will send plots to fgmax_plotdir = \n  ', fgmax_plotdir)
    print('fgreport believes rundir is: ',rundir)

    use_force_dry = False
    if use_force_dry:
        fname_force_dry = os.path.join(input_dir, 'force_dry_init.data')
        print('Using force_dry_init from ', fname_force_dry)

    def savefigp(fname):
        if save_figs:
            fullname = '%s/%s' % (fgmax_plotdir, fname)
            savefig(fullname)
            print('Created ', fullname)
        else:
            print('save_figs = False')


    ### Read in and process the fgmax results from the latest run
    # 
    #
    #xves, yves = -1.2415583333e+02,  4.6981851852e+01 
    #
    #This is the Option 1 at the Ocean Shores Elementary School
    #xves, yves = -124.1549074074074,  46.978796296296295 

    print('outdir = ',outdir)
    t_files = glob.glob(outdir + '/fort.t0*')
    times = []
    for f in t_files:
        lines = open(f,'r').readlines()
        for line in lines:
            if 'time' in line: 
                t = float(line.split()[0])
        times.append(t)
    times.sort()
    print('Output times found: ',times)
    if len(times) > 0:
        t_hours = times[-1] / 3600.
        print('\nfgmax results are presumably from final time: %.1f seconds = %.2f hours'\
             % (times[-1], t_hours))
    else:
        t_hours = nan


    # Read fgmax data:
    fg = fgmax_tools.FGmaxGrid()
    fgmax_input_file_name = outdir + '/fgmax_grids.data'
    print('fgmax input file: \n  %s' % fgmax_input_file_name)
    fg.read_fgmax_grids_data(fgno=fgno, data_file=fgmax_input_file_name)
    fg.read_output(outdir=outdir, indexing='xy')
    fg.tstart_max


    # ### Read B0 from special run:

    if (fgno==1):
        B0_name = 'Hoh'
    elif (fgno==2):
        B0_name = 'Ruby'
    else:
        B0_name = 'Kalaloch'
    fname = input_dir + '/%s_B0.asc' %B0_name
    topoB0 = topotools.Topography()
    topoB0.read(fname, topo_type=3)
    B0 = topoB0.Z
    B0_masked = ma.masked_array(B0, fg.B.mask)
    fg.dz = fg.B - B0_masked
    fg.B0 = B0  #B0_masked

    ###
    print(' ')
    print('min/max dz over fgmax region: %.2f, %.2f' \
            % (fg.dz.min(),fg.dz.max()))
    print('min/max B0 over fgmax region: %.2f, %.2f' \
            % (fg.B0.min(),fg.B0.max()))
    print('min/max B over fgmax region: %.2f, %.2f' \
            % (fg.B.min(),fg.B.max()))
    print(' ')
    ###

    ### TESTING
    ### fg.X and fg.Y are used in the plotting below, so
    ### make sure they are the same as those that come from
    ### reading B0 back in (which was created with Topography Tool)

    if 1:
        print('shape of fg.X and fg.Y were: ',shape(fg.X),shape(fg.Y))
        print('shape of topoB0.X and topoB0.Y were: ',shape(topoB0.X),\
                                                     shape(topoB0.Y))
        print('fg.X max was: ',fg.X.max())
        print('topoB0.X max was: ',topoB0.X.max())
        print('fg.X min was: ',fg.X.min())
        print('topoB0.X min was: ',topoB0.X.min())
        print('fg.Y max was: ',fg.Y.max())
        print('topoB0.Y max was: ',topoB0.Y.max())
        print('fg.Y min was: ',fg.Y.min())
        print('topoB0.Y min was: ',topoB0.Y.min())
        print(' ')


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
        pc = plottools.pcolorcells(fg.X,fg.Y,Z,cmap=cmap,norm=norm)  
        if show_cb:
            cb = colorbar(pc,shrink=0.5,extend='both')
            cb.set_label('meters')
        gca().set_aspect(1./cos(47.74*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20);
    
    figure(figsize=(8,8))
    subplot(211)
    plotZ(fg.B0, show_cb=True)
    title('GeoClaw B0 before quake');

    subplot(212)
    plotZ(fg.B, show_cb=True)
    title('GeoClaw B after quake');
    tight_layout()
    savefigp(B0_name + '_' + 'geoclaw_topo.png')


    onshore = fg.B0 >  0.
    if fg.force_dry_init is not None:
        onshore = logical_or(onshore, fg.force_dry_init)
    offshore = logical_not(onshore)

    fg.h_onshore = ma.masked_where(offshore, fg.h)

    #Use B0 for continuity at the shore
    fg.eta_offshore = ma.masked_where(onshore, fg.B0 + fg.h)

    print('min/max h_onshore over fgmax region: %.2f, %.2f' \
            % (fg.h_onshore.min(),fg.h_onshore.max()))
    print('min/max speed over fgmax region: %.2f, %.2f' \
            % (fg.s.min(),fg.s.max()))
    print(' ')


    #### Plot maximum flow depth

    bounds_depth = array([1e-6,1.0,2.0,4.0,6.0,8.0,10.0,12.0])
    cmap_depth = colors.ListedColormap([[.9,.9,1],[.7,.7,1],[.5,.5,1],\
                  [0,0,1],[1,.7,.7],[1,.4,.4],[1,0,0]])

    # Set color for value exceeding top of range to purple:
    cmap_depth.set_over(color=[1,0,1])

    # Set color for land points without inundation to light green:
    cmap_depth.set_under(color=[.7,1,.7])

    norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)
    
    figure(figsize=(12,8))
    pc = plottools.pcolorcells(fg.X,fg.Y,fg.h_onshore,cmap=cmap_depth,\
         norm=norm_depth)
    cb = colorbar(pc, extend='max', shrink=0.7)
    cb.set_label('meters')
    contour(fg.X, fg.Y, fg.B0, [0], colors='g')
    gca().set_aspect(1./cos(47.74*pi/180.))
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    title('Maximum onshore flow depth over %.2f hours' % t_hours)
    #plot([xves],[yves],'kx',markersize=10)
    savefigp(B0_name + '_' + 'h_onshore.png')

    #### Plot maximum speed
    #bounds_speed = np.array([1e-6,0.5,1.0,1.5,2,3,4.5,6])
    bounds_speed = np.array([1e-6,1,2,4,6,8,10,12])
    cmap_speed = mpl.colors.ListedColormap([[.9,.9,1],[.6,.6,1],\
                 [.3,.3,1],[0,0,1],[1,.8,.8],[1,.6,.6],[1,0,0]])

    # Set color for value exceeding top of range to purple:
    cmap_speed.set_over(color=[1,0,1])

    # Set color for land points without inundation to light green:
    cmap_speed.set_under(color=[.7,1,.7])

    norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)

    figure(figsize=(12,8))
    pc = plottools.pcolorcells(fg.X, fg.Y, fg.s, cmap=cmap_speed,\
         norm=norm_speed)
    cb = colorbar(pc, extend='max', shrink=0.7)
    cb.set_label('m/s')
    contour(fg.X, fg.Y, fg.B0, [0], colors='g')
    gca().set_aspect(1./cos(47.74*pi/180.))
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    title('Maximum speed over %.2f hours' % t_hours)
    #plot([xves],[yves],'kx',markersize=10)
    savefigp(B0_name + '_' + 'speed.png')

    #### Plot maximum eta_offshore
    #bounds_eta = array([0,0.5,1.0,1.5,2,4.0,6.0])
    bounds_eta = array([0,1,2,4,6,10,12])

    cmap_eta = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],
                     [1,.7,.7], [1,.4,.4], [1,0,0]])

    # Set color for value exceeding top of range to purple:
    cmap_eta.set_over(color=[1,0,1])

    norm_eta = colors.BoundaryNorm(bounds_eta, cmap_eta.N)
    
    figure(figsize=(12,8))
    pc = plottools.pcolorcells(fg.X, fg.Y, fg.eta_offshore, cmap=cmap_eta,\
         norm=norm_eta)
    cb = colorbar(pc, extend='max', shrink=0.7)
    cb.set_label('meters')
    contour(fg.X, fg.Y, fg.B0, [0], colors='g')
    gca().set_aspect(1./cos(47.74*pi/180.))
    ticklabel_format(useOffset=False)
    xticks(rotation=20)
    title('Maximum offshore surface eta over %.2f hours' % t_hours)
    #plot([xves],[yves],'kx',markersize=10)
    savefigp(B0_name + '_' + 'eta_offshore.png')

    ### Plots for Google Earth overlays
    # 
    # Tne new version of `kmltools` includes some tools to make png
    # files that display properly on Google Earth.  The png files have
    # no axes and have the dimension and dpi set properly so that there
    # is an integer number of pixels in each grid cell so cell edges are
    # sharp when zooming in. We make three png files and then make a kml
    # file that can be used to open all three.

    if 1:
        kml_dir = fgmax_plotdir + '/' + B0_name +'_' + 'kmlfiles'
        print('Will send kml file and plots to kml_dir = \n  ', kml_dir)
        os.system('mkdir -p %s' % kml_dir);

        #fg.x = fg.X[:,0]
        #fg.y = fg.Y[0,:]

        h_wet_onshore = ma.masked_where(fg.h_onshore==0., fg.h_onshore)
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
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x,fg.y, 
                                    speed,png_filename=png_filename,
                                    dpc=2,cmap=cmap_speed,norm=norm_speed)
        if close_figs: close('all')

        stays_dry = ma.masked_where(fg.h>0., fg.h)
        png_filename = '%s/stays_dry_for_kml.png' % kml_dir
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x,fg.y,
                                    stays_dry, 
                                    png_filename=png_filename,
                                    dpc=2,cmap=cmap_speed,norm=norm_speed)
        if close_figs: close('all')

        ##### Make colorbars for kml files

        kmltools.kml_build_colorbar('%s/colorbar_depth.png' % kml_dir,
                 cmap_depth,norm=norm_depth,label='meters',title='depth',
                 extend='max')
        kmltools.kml_build_colorbar('%s/colorbar_speed.png' % kml_dir,
                 cmap_speed,norm=norm_speed,label='meters / second',
                 title='speed',extend='max')
        kmltools.kml_build_colorbar('%s/colorbar_eta.png' % kml_dir, cmap_eta, 
                 norm=norm_eta,label='meters',title='eta',extend='max')
        if close_figs: close('all')

        ##### Make the kml file to display these three png files
        # Then you can open `fgmax_results_kmlfiles/fgmax_results.kml` 
        # in Google Earth to view them.

        png_files=['h_onshore_max_for_kml.png', 'speed_max_for_kml.png',
                   'stays_dry_for_kml.png',
                   'eta_offshore_max_for_kml.png']
        png_names=['max depth onshore','max speed','stays dry',
                   'eta_offshore']
        cb_files = ['colorbar_depth.png', 'colorbar_speed.png',
                    'colorbar_eta.png']
        cb_names = ['colorbar_depth', 'colorbar_speed',
                    'colorbar_eta']

        name = 'fgmax_%s_%s' % (B0_name,event)
        fname = os.path.join(kml_dir, name+'.kml')
        kmltools.png2kml(png_extent,png_files=png_files,png_names=png_names, 
                         name=name,fname=fname,
                         radio_style=False,
                         cb_files=cb_files,cb_names=cb_names)

        #### Create .kmz file including all plots

        savedir = os.getcwd()
        os.chdir(kml_dir)
        files = glob.glob('*.kml') + glob.glob('*.png')
        print('kmz file will include:')
        for file in files:
            print('    %s' % os.path.split(file)[-1])

            fname_kmz = 'fgmax_results_%s_%s.kmz' % (B0_name,event)
            with zipfile.ZipFile(fname_kmz, 'w') as zip:
                for file in files:
                    zip.write(file) 
                print('Created %s' % os.path.abspath(fname_kmz))
        os.chdir(savedir)

