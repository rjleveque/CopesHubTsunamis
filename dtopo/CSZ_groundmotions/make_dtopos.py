
"""
Convert ground motions into dtopo files
"""

from pylab import *
import os,sys
from clawpack.geoclaw import topotools,kmltools
from scipy.interpolate import LinearNDInterpolator
import obspy
from clawpack.visclaw import animation_tools
from IPython.display import HTML



#CHT = os.environ['CHT'] # path to CopesHubTsunamis directory
CHT = '/Users/rjl/git/CopesHubTsunamis' # or hard-wired
sys.path.insert(0,CHT+'/common_code')  # add to search path
import dtopotools # local version, makes smaller files


# ## Select an event:


event = 'buried-locking-skl16-shallow'
#event = 'buried-random-str10-shallow'
#event = 'buried-locking-str10-deep'
print('Rupture name: ',event)


## Read static displacement

datadir = 'vertical_displacements'
fname_orig = 'vert_displacements_all_xgrid_' + event
path_orig = os.path.join(datadir, fname_orig)
lon,lat,zdisp = loadtxt(path_orig, skiprows=1,usecols=[1,2,3],unpack=True)

ii = argmax(zdisp)
print('max zdisp[%i] = %.1f m at x = %.4f, y = %.4f' \
      % (ii,zdisp[ii],lon[ii],lat[ii]))


# ## Read waveforms

waveforms = obspy.read('time_dependent_zdisp/XGRID_%s_Z.h5' % event)

waveforms = waveforms.sort(['station'])

print('Loaded %i waveforms' % len(waveforms))
tr = waveforms[0]
print('at times ', list(tr.times()))

# ## Create time-dependent dtopo object:


mx = 6*120 + 1  # 30 arcsec
my = 10*120 + 1  # 30 arcsec
x = linspace(-128.5,-122.5,mx)
y  = linspace(40,50,my)
X,Y = meshgrid(x,y)

dtopo = dtopotools.DTopography()
dtopo.X = X
dtopo.Y = Y
dtopo.times = tr.times()
ntimes = len(dtopo.times)
dZshape = (ntimes,X.shape[0],X.shape[1])
dtopo.dZ = empty(dZshape)

points = vstack((lon,lat)).T

for k in range(ntimes):
    # interpolate to uniform grid at each time:
    dispk = [w.data[k] for w in waveforms]
    dz_fcn = LinearNDInterpolator(points, dispk, fill_value=0.)
    dtopo.dZ[k,:,:] = dz_fcn(X,Y)


# ### Save dtopo file for GeoClaw:

fname_dtopo = event + '.dtt3'
dtopo.write(fname_dtopo, 3)
print('Created ',fname_dtopo)

# ### Save dtopo file for instantaneous displacement

dtopo_instant = dtopotools.DTopography()
dtopo_instant.X = X
dtopo_instant.Y = Y
dtopo_instant.times = [0.]
dZshape = (1,X.shape[0],X.shape[1])
dtopo_instant.dZ = empty(dZshape)
dtopo_instant.dZ[0,:,:] = dtopo.dZ[-1,:,:]

fname_dtopo = event + '_instant.dtt3'
dtopo_instant.write(fname_dtopo, 3)
print('Created ',fname_dtopo)


if 1:
    # ## Make plots and animations:


    coast = load('/Users/rjl/git/clawpack/geoclaw/scratch/pacific_shorelines_east_4min.npy')
    x_coast = coast[:,0] - 360.
    y_coast = coast[:,1]


    dZ_max = abs(dtopo.dZ).max()
    cmax_dZ = round(dZ_max)
    print('maximum |dZ| = %.2f m, setting cmax_dZ = %.1f' % (dZ_max,cmax_dZ))


    # plot at final time:
    fig,ax = subplots(figsize=(6,7))
    plot(x_coast,y_coast,'g')
    axis([-128.5,-122,40,50])
    t = dtopo.times.max()
    dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
    title('Final vertical displacement\n%s' % event)
    grid(True);
    fname = '%s_final.png' % event
    savefig(fname)
    print('Created ',fname)


    # make animation

    figs = []
    for k,t in enumerate(dtopo.times):
        fig,ax = plt.subplots(figsize=(6,7))
        plot(x_coast,y_coast,'g')
        axis([-128.5,-122,40,50])
        dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
        title('%i seconds\n%s' % (t,event))
        grid(True)
        figs.append(fig)
        close(fig)


    anim = animation_tools.animate_figs(figs, figsize=(7,8))

    # make mp4 file:
    fname = event + '.mp4'
    animation_tools.make_mp4(anim, file_name=fname)
    print('Created ',fname)


if 1:
    # ## plot transect of dz at a sequence of times:

    y0 = 47.5
    j = where(y<y0)[0].max()
    for k in range(6,15):
        plot(x,dtopo.dZ[k,j,:],label='%.0fs' % dtopo.times[k])
    legend(loc='upper left')
    title('Transect at y = %.1f\n%s' % (y0,event))
    grid(True)
    fname = event + '_y%s.png' % str(y0).replace('.','-')
    savefig(fname)
    print('Created ',fname)
