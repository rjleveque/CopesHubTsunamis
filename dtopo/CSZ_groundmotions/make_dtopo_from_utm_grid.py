
from pylab import *
import os,sys
from clawpack.geoclaw import topotools,kmltools
from scipy.interpolate import LinearNDInterpolator

#CHT = os.environ['CHT'] # path to CopesHubTsunamis directory
CHT = '/Users/rjl/git/CopesHubTsunamis' # or hard-wired
sys.path.insert(0,CHT+'/common_code')  # add to search path
import dtopotools # local version, makes smaller files

#defdir = 'vertical_displacements' # first 18 events
#defdir = 'vertical_displacements_sensitivity'  # 3 new events 6/20/24
defdir = 'vertical_displacements_FrontalThrust'  # new events 11/7/24

#fname_orig = 'vert_displacements_all_xgrid_buried-random-str10-shallow'
fname_orig = 'vert_displacements_all_xgrid_ft-random-mur13-deep'
#fname_orig = 'vert_displacements_all_xgrid_buried-random-mur13-deep-nosub'
#fname_orig = 'vert_displacements_all_xgrid_buried-random-mur13-deep-nosub-homogeneous-halfspace-noQ'

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
dtopo.times = array([0])
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

# simplify name for sensitivity study:
fname_orig = fname_orig.replace('vert_displacements_all_xgrid_','')
fname_orig = fname_orig + '_instant'  # since these are static displacements

fname = fname_orig + '.dtt3'
path_new = os.path.join(defdir, fname)
dtopo.write(path_new, 3)
print('Created ', path_new)


if 1:
    # Make plots and animations:


    coast = load('/Users/rjl/git/clawpack/geoclaw/scratch/pacific_shorelines_east_4min.npy')
    x_coast = coast[:,0] - 360.
    y_coast = coast[:,1]


    dZ_max = abs(dtopo.dZ).max()
    #cmax_dZ = round(dZ_max)
    cmax_dZ = 4. # to compare to Audrey's plot
    print('maximum |dZ| = %.2f m, setting cmax_dZ = %.1f' % (dZ_max,cmax_dZ))

    # plot at final time (only time for instant):
    fig,ax = subplots(figsize=(6,7))
    plot(x_coast,y_coast,'g')
    axis([-128.5,-122,40,50])
    t = dtopo.times.max()
    dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
    title('Final vertical displacement\n%s' % fname_orig)
    grid(True);
    fname = '%s/%s.png' % (defdir,fname_orig)
    savefig(fname)
    print('Created ',fname)


