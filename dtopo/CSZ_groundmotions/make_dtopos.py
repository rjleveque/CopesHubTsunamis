
"""
Convert ground motions into dtopo files
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools
from scipy.interpolate import LinearNDInterpolator
import obspy
from clawpack.visclaw import animation_tools
from IPython.display import HTML



#CHT = os.environ['CHT'] # path to CopesHubTsunamis directory
CHT = '/Users/rjl/git/CopesHubTsunamis' # or hard-wired
sys.path.insert(0,CHT+'/common_code')  # add to search path
import dtopotools # local version, makes smaller files

zdisp_dir = 'time_dependent_zdisp_FrontalThrust'
dtopo_dir = 'dtopofiles_FrontalThrust'
os.system('mkdir -p %s' % dtopo_dir)
print('Output will go in %s' % dtopo_dir)

# Select an event:

if 0:
    #event = 'buried-locking-str10-deep'
    event = 'ft-locking-mur13-deep'
    events = [event]

if 1:
    # all events in zdisp_dir:
    #files = glob.glob('vertical_displacements_FrontalThurst/*')
    #files = glob.glob('time_dependent_zdisp/*')
    files = glob.glob('%s/*' % zdisp_dir)
    #files = ['time_dependent_zdisp_FrontalThrust/XGRID_ft-locking-mur13-deep_Z.h5'])
    events = []
    for f in files:
        f = os.path.split(f)[-1]
        #events.append(f.replace('vert_displacements_all_xgrid_', ''))
        f = f.replace('XGRID_', '')
        f = f.replace('_Z.h5', '')
        #if 'locking-mur13' in f:
        if 1:
            events.append(f)
    print(events)


for event in events:
    print('Rupture name: ',event)

    if 0:
        # Read static displacement  (not needed)
        datadir = 'vertical_displacements' # might be wrong
        fname_orig = 'vert_displacements_all_xgrid_' + event
        path_orig = os.path.join(datadir, fname_orig)
        lon,lat,zdisp = loadtxt(path_orig, skiprows=1,usecols=[1,2,3],unpack=True)
        points = vstack((lon,lat)).T

        ii = argmax(zdisp)
        print('max zdisp[%i] = %.1f m at x = %.4f, y = %.4f' \
              % (ii,zdisp[ii],lon[ii],lat[ii]))

    points = loadtxt('lonlat_points.txt')  # created by get_lonlat.py

    # Read waveforms

    waveforms = obspy.read('%s/XGRID_%s_Z.h5' % (zdisp_dir,event))

    waveforms = waveforms.sort(['station'])

    print('Loaded %i waveforms' % len(waveforms))
    tr = waveforms[0]
    print('at times ', list(tr.times()))

    # Create time-dependent dtopo object:


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


    for k in range(ntimes):
        # interpolate to uniform grid at each time:
        dispk = [w.data[k] for w in waveforms]
        dz_fcn = LinearNDInterpolator(points, dispk, fill_value=0.)
        dtopo.dZ[k,:,:] = dz_fcn(X,Y)

    if 0:
        # truncate to shorter time period if motion stops early
        # diff_to_end didn't go below 0.34 for one test, so dropping this...
        diff_to_end = zeros(ntimes)
        for k in range(ntimes):
            diff_to_end[k] = abs(dtopo.dZ[k,:,:] - dtopo.dZ[-1,:,:]).max()

        print('diff_to_end = ',diff_to_end)
        kend = ntimes
        for k in range(kend-1,0,-1):
            if diff_to_end[k:].max() < 0.001:
                kend = k
        print('original ntimes = %i,  kend = %i' % (ntimes,kend))
        dtopo.dZ = dtopo.dZ[:kend,:,:]
        dtopo.times = dtopo.times[:kend]
        ntimes = len(dtopo.times)
        print('new ntimes = %i' % ntimes)

    # Save dtopo file for GeoClaw:

    fname_dtopo = '%s/%s.dtt3' % (dtopo_dir,event)
    dtopo.write(fname_dtopo, 3)
    print('Created ',fname_dtopo)

    # Save dtopo file for instantaneous displacement

    dtopo_instant = dtopotools.DTopography()
    dtopo_instant.X = X
    dtopo_instant.Y = Y
    dtopo_instant.times = [0.]
    dZshape = (1,X.shape[0],X.shape[1])
    dtopo_instant.dZ = empty(dZshape)
    dtopo_instant.dZ[0,:,:] = dtopo.dZ[-1,:,:]

    fname_dtopo = '%s/%s_instant.dtt3' % (dtopo_dir,event)
    dtopo_instant.write(fname_dtopo, 3)
    print('Created ',fname_dtopo)


    if 1:
        # Make plots and animations:


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
        fname = '%s/%s_final.png' % (dtopo_dir,event)
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
        fname = '%s/%s.mp4' % (dtopo_dir,event)
        animation_tools.make_mp4(anim, file_name=fname)
        print('Created ',fname)


    if 1:
        for y0 in [44,45,46,47,48]:
            # plot transect of dz at a sequence of times:

            time_interval = 4
            time_indices = range(time_interval, len(dtopo.times), time_interval)
            #times = dtopo.times[interval::interval]
            figure(figsize=(10,7))
            #y0 = 47.5
            j = where(y<y0)[0].max()
            for k in time_indices:
                plot(x,dtopo.dZ[k,j,:],label='%.0fs' % dtopo.times[k])
            legend(loc='upper left')
            title('Transect at y = %.1f\n%s' % (y0,event))
            grid(True)
            fname = '%s/%s_y%s.png' % (dtopo_dir,event, str(y0).replace('.','-'))
            savefig(fname)
            print('Created ',fname)
