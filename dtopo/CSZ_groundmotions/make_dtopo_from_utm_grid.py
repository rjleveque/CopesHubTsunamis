
from pylab import *
import os,sys
from clawpack.geoclaw import topotools,dtopotools,kmltools
from scipy.interpolate import LinearNDInterpolator

defdir = 'vertical_displacements'
#fname_orig = 'vert_displacements_all_xgrid_buried-random-str10-shallow'
fname_orig = 'vert_displacements_all_xgrid_buried-locking-skl16-shallow'
path_orig = os.path.join(defdir, fname_orig)
lon,lat,z = loadtxt(path_orig, skiprows=1,usecols=[1,2,3],unpack=True)

points = vstack((lon,lat)).T
z_fcn = LinearNDInterpolator(points, z, fill_value=0.)

mx = 6*120 + 1  # 30 arcsec
my = 10*120 + 1  # 30 arcsec
x = linspace(-128.5,-122.5,mx)
y  = linspace(40,50,my)
X,Y = meshgrid(x,y)
Z = z_fcn(X,Y)

dtopo = dtopotools.DTopography()
dtopo.dZ = Z
dtopo.X = X
dtopo.Y = Y
dtopo.times = [0]
ntimes = len(dtopo.times)
#dtopo.delta = (dx,dy)
dZshape = (ntimes,X.shape[0],X.shape[1])
dZ = empty(dZshape)

if ntimes==1:
    dZ[0,:,:] = Z
else:
    for k in range(ntimes):
        dZ[k,:,:] = k/(ntimes-1) * Z

dtopo.dZ = dZ

fname = fname_orig + '.dtt3'
path_new = os.path.join(defdir, fname)
dtopo.write(fname, 3)
print('Created ', fname)


