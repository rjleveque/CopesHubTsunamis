from pylab import *
from clawpack.geoclaw import dtopotools
import contextlib, os, glob

if 0:
    datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/jey_250614/'
    lon,lat,z = loadtxt(datadir+'poly3d_uz_resps_Poly3D_2012_finer_spacing.xyz',
                    unpack=True)

if 0:
    datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/jey_250616/output_files_shared_06162025/'
    fname = 'poly3d_uz_resps_source_saved_locking_mur13_deep.out'
    fname_dtopo = fname.replace('poly3d_uz_resps_source_saved_','')
    fname_dtopo = datadir + fname.replace('.out','_poly3d_5km.dtt3')
    lon,lat,z = loadtxt(datadir+fname, unpack=True)

if 0:
    datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/audrey_250616/input_vertical_displacements/vertical_disp_poly3d_5x5km_resolution/'
    fname = 'buried_random_mur13_deep_Poly3D_vertdisp.txt'
    fname_dtopo = datadir + fname.replace('.txt','_5km.dtt3')
    lon,lat,z = loadtxt(datadir+fname, unpack=True, skiprows=1)


datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/jey_250616/output_files_shared_06162025/'

dtopodir = 'dtopofiles_poly3d_5km/'
os.system('mkdir -p %s' % dtopodir)

with contextlib.chdir(datadir):
    fnames = glob.glob('poly3d_uz_resps_source_saved*.out')

for fname in fnames:
    fname_dtopo = fname.replace('poly3d_uz_resps_source_saved_','buried_')
    fname_dtopo = dtopodir + fname_dtopo.replace('.out','_poly3d_5km_instant.dtt3')
    lon,lat,z = loadtxt(datadir+fname, unpack=True)

    mx = 121
    my = 201
    dz = reshape(z, (mx,my), order='F').T

    x = linspace(lon.min(),lon.max(),mx)
    y = linspace(lat.min(),lat.max(),my)
    X,Y = meshgrid(x,y,indexing='xy')

    dtopo = dtopotools.DTopography()
    dtopo.X = X
    dtopo.Y = Y
    dtopo.dZ = array([zeros(dz.shape), flipud(dz)])
    dtopo.times = array([0,1])

    dtopo.write(fname_dtopo, dtopo_type=3)
    print('Created ',fname_dtopo)
