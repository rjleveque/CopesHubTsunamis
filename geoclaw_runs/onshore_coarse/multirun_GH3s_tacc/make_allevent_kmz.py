
"""
Process fgmax grid results and make kmz file showing max depth for all
36 events with single color bar.

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


try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Set CHT enviornment variable to repository top")

# location for big files for different computer environments:
this_dir = os.getcwd()
HOME = os.environ['HOME']

if 'rjl/git' in this_dir:
    computer = 'rjl-laptop'
    scratch_dir = this_dir.replace('rjl/git', 'rjl/scratch')

elif '/mmfs1/home' in this_dir:
    computer = 'hyak'
    scratch_dir = this_dir.replace('/mmfs1/home', '/gscratch/tsunami')

elif '/home1' in this_dir:
    computer = 'tacc'
    #scratch_dir = this_dir.replace('/home1', '/scratch')
    try:
        SCRATCH = os.environ['SCRATCH']
        scratch_dir = this_dir.replace(HOME, SCRATCH)
    except:
        scratch_dir = this_dir  # if $SCRATCH not set

else:
    computer = 'unknown'
    scratch_dir = this_dir

use_force_dry = False
if use_force_dry:
    fname_force_dry = os.path.join(input_dir, 'force_dry_init.data')
    print('Using force_dry_init from ', fname_force_dry)


def load_fgmax(outdir,fgno,fname_B0=None):
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

    if fname_B0:
        topoB0 = topotools.Topography()
        topoB0.read(fname_B0, topo_type=3)
        B0 = topoB0.Z
        B0_masked = ma.masked_array(B0, fg.B.mask)
        fg.B0 = B0_masked
    else:
        fg.B0 = fg.B
        print('Assuming B0 = fg.B')

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



def make_kmz_plots(fgno, fg, plotdir, event):
    # ## Plots for Google Earth overlays

    kml_dir = plotdir + '/kmlfiles'
    print('Will send kml file and plots to kml_dir = \n  ', kml_dir)
    os.system('mkdir -p %s' % kml_dir);


    if 1: # usually
        othercondition = False    #onshore determined only by fg.B0 > 0

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


    h_wet_onshore = ma.masked_where(fg.h_onshore==0., fg.zeta_onshore)
    print('fg.x, fg.y shapes: ',fg.x.shape, fg.y.shape)
    print('+++ h_wet_onshore.shape = ',h_wet_onshore.shape)
    png_filename = kml_dir+'/%s_h_onshore_max_for_kml.png' % event
    fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                     h_wet_onshore,
                                                     png_filename=png_filename,
                                                     dpc=2, cmap=cmap_depth, norm=norm_depth)
    if close_figs: close('all')

    if 0:

        png_filename=kml_dir+'/%s_eta_offshore_max_for_kml.png' % event
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                         fg.eta_offshore,
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_eta, norm=norm_eta)
        if close_figs: close('all')



        speed = ma.masked_where(fg.h==0., fg.s)
        png_filename = '%s/%s_speed_max_for_kml.png' % (kml_dir,event)
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y, speed,
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_speed, norm=norm_speed)
        if close_figs: close('all')


        stays_dry = ma.masked_where(fg.h>0., fg.h)
        png_filename = '%s/%s_stays_dry_for_kml.png' % (kml_dir,event)
        fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(fg.x, fg.y,
                                                         stays_dry,
                                                         png_filename=png_filename,
                                                         dpc=2, cmap=cmap_speed, norm=norm_speed)
        if close_figs: close('all')

    return png_filename, png_extent

def make_colorbar_kmls(plotdir, name_kmz, png_extent):
    # ### Make colorbars for kml files

    kml_dir = plotdir + '/kmlfiles'

    if 1:
        kmltools.kml_build_colorbar('%s/colorbar_depth.png' % kml_dir, cmap_depth,
                                   norm=norm_depth, label='meters', title='depth', extend='max')

        cb_files = ['colorbar_depth.png']
        cb_names = ['colorbar_depth']

    if 0:
        kmltools.kml_build_colorbar('%s/colorbar_speed.png' % kml_dir, cmap_speed,
                                   norm=norm_speed, label='meters / second', title='speed', extend='max')
        kmltools.kml_build_colorbar('%s/colorbar_eta.png' % kml_dir, cmap_eta,
                                   norm=norm_eta, label='meters', title='eta', extend='max')
        cb_files = ['colorbar_depth.png', 'colorbar_speed.png',
                    'colorbar_eta.png']
        cb_names = ['colorbar_depth', 'colorbar_speed',
                    'colorbar_eta']

    if close_figs: close('all')
    return cb_files, cb_names



def make_main_kml(name_kmz, png_extent, png_files, cb_files, cb_names):


    png_names=[ss for s in png_files]


    # make kml file including everything (this will be only kml in kmz file):

    name = name_kmz
    fname = os.path.join(kml_dir, name+'.kml')

    kmltools.png2kml(png_extent, png_files=png_files, png_names=png_names,
                     name=name, fname=fname,
                     radio_style=False,
                     cb_files=cb_files, cb_names=cb_names)


def make_all_kmz_plots(events, outdirs, plotdir, name_kmz):

    kml_dir = plotdir + '/kmlfiles'
    os.system('mkdir -p %s' % kml_dir)

    path_kmz = os.path.join(plotdir, name_kmz + '.kmz')

    fgno = 1

    png_files = []
    png_names = []
    for event,outdir in zip(events,outdirs):

        fname_B0 = None

        print('Will read fgmax results from outdir = \n  ', outdir)
        fg, t_hours = load_fgmax(outdir,fgno,fname_B0)

        png_filename, png_extent = make_kmz_plots(fgno, fg, plotdir, event)

        png_filename = os.path.split(png_filename)[-1]
        #print('+++ png_filename = ',png_filename)
        png_files.append(png_filename)
        png_names.append(event)
        if close_figs: close('all')

    # make colorbar(s):
    cb_files, cb_names = make_colorbar_kmls(plotdir, name_kmz, png_extent)

    # make main kml file (will be only kml in kmz file, along with pngs):

    name = name_kmz
    fname = os.path.join(kml_dir, name+'.kml')

    kmltools.png2kml(png_extent, png_files=png_files, png_names=png_names,
                     name=name, fname=fname,
                     radio_style=False,
                     cb_files=cb_files, cb_names=cb_names)

    if close_figs: close('all')

    # make kmz file including everything:

    savedir = os.getcwd()
    os.chdir(kml_dir)
    files = glob.glob('*.kml') + glob.glob('*.png')
    print('kmz file will include:')
    for file in files:
        print('    %s' % os.path.split(file)[-1])

    with zipfile.ZipFile(path_kmz, 'w') as zipf:
        for file in files:
            zipf.write(file)

    #path_kmz = os.path.join(fgmax_plotdir, name_kmz)
    #shutil.move(name_kmz, path_kmz)
    print('Created %s' % os.path.abspath(path_kmz))
    if 0:
        print('kml and png files are in: %s' % kml_dir)
    else:
        shutil.rmtree(kml_dir)
    os.chdir(savedir)


# List of all events from CoPes Hub ground motions:
# For naming and numbering convention, see
#   https://depts.washington.edu/ptha/CHTuser/dtopo/groundmotions/
# Note that there may be different versions of these events, so which
# version is used depends on what directory dtopo_dir points to

try:
    import CHTtools
    all_events = CHTtools.all_events()  # returns the list of 36 events
except:
    # if CHTtools isn't found, construct the list of 36 events here:

    depths = ['D','M','S']

    # buried_locking events:
    all_events = [f'BL10{depth}' for depth in depths] \
               + [f'BL13{depth}' for depth in depths] \
               + [f'BL16{depth}' for depth in depths] \

    all_events += [e.replace('L','R') for e in all_events]  # add random events
    all_events += [e.replace('B','F') for e in all_events]  # add ft events

    all_events.sort()


# if dtopo_dir points to a directory that has instantaneous versions
# (static displacement rather than kinematic time-dependent)
# then you could set `instant = True` to use these, provided they
# have file names such as BL10D_instant.dtt3 (with the same numbering 1-36):

instant = False

if instant:
    all_events = [e+'_instant' for e in all_events]

events = all_events[:4]

name_kmz = 'GH3s_4events'

if __name__== '__main__':

    # location for outdirs:
    this_dir = os.getcwd()

    # Randy's laptop:
    scratch_dir = this_dir.replace('git/CopesHubTsunamis/geoclaw_runs', \
                                   'scratch/CHT_runs')
    #scratch_dir = '/Users/rjl/tests/CHT_runs'

    # for hyak:
    scratch_dir = scratch_dir.replace('/mmfs1/home', '/gscratch/tsunami')

    runs_dir = os.path.abspath(scratch_dir)

    runs_dir = os.path.abspath('hyak_geoclaw_outputs')  # on laptop

    print('+++ this_dir = ',this_dir)
    print('+++ runs_dir = ',runs_dir)

    all_models = []

    if 1:
        all_models = all_models + \
            ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
             'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']
        name_kmz = 'coarse_hmax_GH3s_buried'

    if 0:
        all_models = all_models + \
            ['ft-locking-mur13', 'ft-locking-skl16', 'ft-locking-str10',
             'ft-random-mur13',  'ft-random-skl16',  'ft-random-str10']
        name_kmz = 'coarse_hmax_GH3s_ft'

    if len(all_models) == 12:
        # including both buried and ft:
        name_kmz = 'coarse_hmax_GH3s'

    models = all_models
    #models = all_models[:3]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

    events.sort()

    #events = events[:3]

    instant = False
    if instant:
        events = [e+'_instant' for e in events]

    if 0:
        #events = ['ft-locking-mur13-deep']
        events = ['buried-locking-mur13-deep']


    outdirs = ['%s/geoclaw_outputs/_output_%s' % (runs_dir, event) \
                for event in events]

    plotdir = '%s/geoclaw_plots' % runs_dir

    make_all_kmz_plots(events, outdirs, plotdir, name_kmz)
