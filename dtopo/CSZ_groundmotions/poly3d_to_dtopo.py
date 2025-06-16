from pylab import *
from clawpack.geoclaw import dtopotools

datadir = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/jey_250614/'
lon2,lat2,z2 = loadtxt(datadir+'poly3d_uz_resps_Poly3D_2012_finer_spacing.xyz',
                unpack=True)

mx = 121
my = 201
dz = reshape(z2, (mx,my), order='F').T

x = linspace(lon2.min(),lon2.max(),mx)
y = linspace(lat2.min(),lat2.max(),my)
X,Y = meshgrid(x,y,indexing='xy')

dtopo = dtopotools.DTopography()
dtopo.X = X
dtopo.Y = Y
dtopo.dZ = array([zeros(dz.shape), flipud(dz)])
dtopo.times = array([0,1])

fname = datadir + 'poly3d_uz_resps_Poly3D_2012_finer_spacing.dtt3'
dtopo.write(fname, dtopo_type=3)
