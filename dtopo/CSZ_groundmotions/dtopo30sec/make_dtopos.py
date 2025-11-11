
"""
Convert ground motions into dtopo files

15-arcsecond resolution with new naming convention, e.g.
     FL13D = ft-locking-mur13-deep
dt = 10 seconds
format '%.2f' (cm resolution)
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools,dtopotools
from scipy.interpolate import LinearNDInterpolator
import obspy
from clawpack.visclaw import animation_tools
from IPython.display import HTML
#from clawpack.clawutil import util


#CHT = os.environ['CHT'] # path to CopesHubTsunamis directory
#CHT = '/Users/rjl/git/CopesHubTsunamis' # or hard-wired
#sys.path.insert(0,CHT+'/common_code')  # add to search path
#import dtopotools # local version, makes smaller files


dtopo_dir = './dtopofiles'
os.system('mkdir -p %s' % dtopo_dir)
print('Output will go in %s' % dtopo_dir)

def shortname(event):
    """
    convert original event name like ft-locking-mur13-deep
    to new short name like FL13D
    """

    if len(event) == 5:
        # assume it is already a short name:
        return event

    tokens = event.replace('_','-').split('-')
    newname = ''
    if tokens[0] == 'buried':
        newname += 'B'
    elif tokens[0] == 'ft':
        newname += 'F'
    if tokens[1] == 'locking':
        newname += 'L'
    elif tokens[1] == 'random':
        newname += 'R'
    newname += tokens[2][-2:]  # year
    if tokens[3] == 'deep':
        newname += 'D'
    elif tokens[3] == 'middle':
        newname += 'M'
    elif tokens[3] == 'shallow':
        newname += 'S'

    if len(newname) != 5:
        print(f'*** problem converting {event}, newname = {newname}')

    return newname

# Select an event:

if 0:
    event = 'buried-random-mur13-deep'
    #event = 'ft-locking-mur13-deep'
    events = [event]

if 1:
    depths = ('deep','middle','shallow')
    events = [f'buried-locking-str10-{depth}' for depth in depths]
    events += [f'buried-locking-mur13-{depth}' for depth in depths]
    events += [f'buried-locking-skl16-{depth}' for depth in depths]
    events += [e.replace('buried','ft') for e in events]
    events += [e.replace('locking','random') for e in events]

if 0:
    # all events in zdisp_dir:
    zdisp_dir = '../time_dependent_zdisp'  # for buried ruptures
    files = glob.glob('%s/*' % zdisp_dir)
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

#events = events[:1]
events_short = [shortname(event) for event in events]
print('Using these events: ',events_short)

for event in events:
    event_short = shortname(event)
    print(f'Rupture name: {event} or {event_short}')

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

    points = loadtxt('../lonlat_points.txt')  # created by get_lonlat.py

    # Read waveforms
    if event[:2] == 'bu':
        zdisp_dir = '../time_dependent_zdisp'  # for buried ruptures
    elif event[:2] == 'ft':
        zdisp_dir = '../time_dependent_zdisp_FrontalThrust'
    waveforms = obspy.read('%s/XGRID_%s_Z.h5' % (zdisp_dir,event))

    waveforms = waveforms.sort(['station'])

    print('Loaded %i waveforms' % len(waveforms))
    tr = waveforms[0]
    print('at times ', array(tr.times()))

    # Create time-dependent dtopo object:


    if 0:
        mx = 6*240 + 1  # 15 arcsec
        my = 10*240 + 1  # 15 arcsec
        #mx = 6*120 + 1  # 30 arcsec
        #my = 10*120 + 1  # 30 arcsec
        x = linspace(-128.5,-122.5,mx)
        y  = linspace(40,50,my)

    # increase x range a bit, and align with 1/15" cell centers
    #mx = int(6.8*240 + 1)  # 15 arcsec
    #my = 10*240 + 1  # 15 arcsec
    mx = 6*120 + 1  # 30 arcsec
    my = 10*120 + 1  # 30 arcsec
    dx2 = 0.5 * 15/3600.
    x = linspace(-128.8+dx2,-122.+dx2,mx)
    y  = linspace(40+dx2,50+dx2,my)
    X,Y = meshgrid(x,y)

    dtopo = dtopotools.DTopography()
    dtopo.x = x
    dtopo.y = y
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

    fname_dtopo = '%s/%s.dtt3' % (dtopo_dir,event_short)
    dtopo.write(fname_dtopo, 3, dZ_format='%.2f')  # cm precision
    print('Created ',fname_dtopo)

    # Save dtopo file for instantaneous displacement

    dtopo_instant = dtopotools.DTopography()
    dtopo_instant.X = X
    dtopo_instant.Y = Y
    dtopo_instant.times = [1.]
    dZshape = (1,X.shape[0],X.shape[1])
    dtopo_instant.dZ = empty(dZshape)
    dtopo_instant.dZ[0,:,:] = dtopo.dZ[-1,:,:]

    fname_dtopo = '%s/%s_instant.dtt3' % (dtopo_dir,event_short)
    dtopo_instant.write(fname_dtopo, 3)
    print('Created ',fname_dtopo)


    if 0:
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
        title('Final vertical displacement\n%s' % event_short)
        grid(True);
        fname = '%s/%s_final.png' % (dtopo_dir,event_short)
        savefig(fname)
        print('Created ',fname)


        # make animation

        figs = []
        for k,t in enumerate(dtopo.times):
            fig,ax = plt.subplots(figsize=(6,7))
            plot(x_coast,y_coast,'g')
            axis([-128.5,-122,40,50])
            dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
            title('%i seconds\n%s' % (t,event_short))
            grid(True)
            figs.append(fig)
            close(fig)


        anim = animation_tools.animate_figs(figs, figsize=(7,8))

        # make mp4 file:
        fname = '%s/%s.mp4' % (dtopo_dir,event_short)
        animation_tools.make_mp4(anim, file_name=fname)
        print('Created ',fname)


    if 0:
        # make transect plots
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
            title('Transect at y = %.1f\n%s' % (y0,event_short))
            grid(True)
            fname = '%s/%s_y%s.png' % (dtopo_dir,event_short,
                                       str(y0).replace('.','-'))
            savefig(fname)
            print('Created ',fname)
