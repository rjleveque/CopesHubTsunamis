from pylab import *
from clawpack.geoclaw import dtopotools
import glob,os
import contextlib

datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/audrey_250616/input_vertical_displacements/vertical_disp_3dsim_5x5km_resolution/'

#fname = 'buried_locking_mur13_deep_3D_vertdisp.txt'


with contextlib.chdir(datadir):
    fnames = glob.glob('*.txt')

#print(fnames)

for fname in fnames:

    lon = []
    lat = []
    z = []
    lines = open('%s/%s' % (datadir,fname)).readlines()
    for line in lines[1:]:
        tokens = line.split()
        lon.append(float(tokens[0]))
        lat.append(float(tokens[1]))
        try:
            zj = float(tokens[2])
        except:
            zj = 0.
        z.append(zj)

    lon = asarray(lon)
    lat = asarray(lat)
    z = asarray(z)

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

    fname = datadir + fname.replace('_vertdisp.txt','_5km_instant.dtt3')
    dtopo.write(fname, dtopo_type=3)
    print('Created ',fname)
