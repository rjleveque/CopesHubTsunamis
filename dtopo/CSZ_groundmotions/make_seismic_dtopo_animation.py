
"""
Make animation of original seismic ground motions together with
KinOkada tsunami generation.
"""

from pylab import *
import os,sys,glob
from scipy.interpolate import LinearNDInterpolator
import obspy
from clawpack.geoclaw import topotools,kmltools,dtopotools,fgout_tools
from clawpack.visclaw import animation_tools, colormaps, geoplot, plottools
import OkadaFromSubsampledSlipNosub


from clawpack.clawutil.util import fullpath_import
CHTuser = os.environ['CHTuser']
CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')

CHT = os.environ['CHT']

#zdisp_dir = 'time_dependent_zdisp_FrontalThrust'
zdisp_dir = 'time_dependent_zdisp'
dtopo_dir = 'dtopo30sec_nosubevents_kinokada/dtopofiles'
plot_dir = 'dtopo30sec_nosubevents_kinokada/plots'
os.system(f'mkdir -p {plot_dir}')

if 1:
    #event = 'buried-locking-str10-deep'
    event = 'buried-random-mur13-deep'  # vZ hardwired below for this event
    events = [event]

if 0:
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

    # Read time-dependent dtopo object:
    event_kinokada = CHTtools.shortname(event)
    print(f'KinOkada event: {event_kinokada}')

    print(f'Loading {event_kinokada}...')
    # just use megathrust even for ft events
    fault = OkadaFromSubsampledSlipNosub.load_fault('M',event_kinokada)
    fault.rupture_type = 'dynamic'

    dtopo = dtopotools.DTopography(f'{dtopo_dir}/{event_kinokada}.dtt3',3)

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

    #waveforms = obspy.read('%s/XGRID_%s_Z.h5' % (zdisp_dir,event))

    # velocities
    waveforms = obspy.read('XGRID_buried-random-mur13-deepvel_samp5s_Z.h5')

    waveforms = waveforms.sort(['station'])

    print('Loaded %i waveforms' % len(waveforms))
    tr = waveforms[0]
    waveform_times = array([float(t) for t in tr.times()])
    print('at times ', waveform_times)
    ntimes = len(waveform_times)

    print('dtopo.times = ',dtopo.times)

    # need to difference waveforms in time to estimate velocities,
    # but need finer temporal resolution to do so.
    # For now, try setting up structure with dZ rather than vZ:

    X = dtopo.X
    Y = dtopo.Y
    #seismic_dZ = zeros((ntimes,X.shape[0],X.shape[1]))
    seismic_vZ = zeros((ntimes,X.shape[0],X.shape[1]))

    for k in range(ntimes):
        # interpolate to uniform grid at each time:
        dispk = [w.data[k] for w in waveforms]
        vz_fcn = LinearNDInterpolator(points, dispk, fill_value=0.)
        seismic_vZ[k,:,:] = vz_fcn(X,Y)

    print(f'Created seismic_vZ with shape {seismic_vZ.shape}')



    #coast = load('/Users/rjl/git/clawpack/geoclaw/scratch/pacific_shorelines_east_4min.npy')
    coast = load('CSZ_coast.npy')
    x_coast = coast[:,0]
    y_coast = coast[:,1]


    dZ_max = abs(dtopo.dZ).max()
    cmax_dZ = round(dZ_max)
    print('maximum |dZ| = %.2f m, setting cmax_dZ = %.1f' % (dZ_max,cmax_dZ))

    if 0:

        # plot at final time:
        fig,ax = subplots(figsize=(6,7))
        plot(x_coast,y_coast,'g')
        axis('scaled')
        axis([-128.5,-122,40,50])
        t = dtopo.times.max()
        dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
        title(f'Final vertical displacement\n{event_kinokada} KinOkada')
        grid(True);
        fname = '%s/%s_final.png' % (plot_dir,event_kinokada)
        savefig(fname)
        print('Created ',fname)



    #cmap_slip = colormaps.white_red
    cmap_slip = colormaps.make_colormap({0:[0,0,0,0], 1:[1,1,0,1]})

    #cmap_v = colormaps.white_blue
    #cmap_v = colormaps.make_colormap({0: [0,0,0,0], 1:[0,0,1,0.5]})
    cmap_v = colormaps.make_colormap({0: [0,0,0,0], 1:'purple'})
    clim_v = 0.05*seismic_vZ.max()

    if 1:

        # make animation of slip/vZ (left) with dZ (right)

        cmap_dZ = colormaps.blue_white_red
        cmax_dZ = 10

        figs = []
        for k,t in enumerate(dtopo.times):
            kv = 2*k
            assert waveform_times[kv] == t, \
                    f'*** waveform_time {waveform_times[kv]} != {t} = t'
            fig,axs = plt.subplots(1,2,figsize=(10,8))
            for kax,ax in enumerate(axs):
                ax.plot(x_coast,y_coast,'g',linewidth=0.7)
                if kax==0:
                    #pc = ax.pcolormesh(X,Y,seismic_dZ[k,:,:],cmap=cmap)
                    #pc.set_clim(-cmax_dZ,cmax_dZ)
                    #ax.set_title(f'dZ at {t} seconds\n{event}\nSeismic, with subevents')
                    fault.plot_subfaults(axes=ax,slip_color=True,
                                         cmap_slip=cmap_slip,
                                         slip_time=t,plot_box=False,cmin_slip=0.01,
                                         colorbar_shrink=0.7)
                    pc = ax.pcolormesh(X,Y,seismic_vZ[kv,:,:],cmap=cmap_v)
                    pc.set_clim(0,clim_v)
                    ax.set_title(f'Velocity vZ at {t} seconds\n{event}\nWith KinOkada slip')
                if kax==1:
                    #dtopo.plot_dZ_colors(t,axes=ax,cmax_dZ=cmax_dZ,dZ_interval=100)
                    pc = ax.pcolormesh(dtopo.X,dtopo.Y,dtopo.dZ[k,:,:],cmap=cmap_dZ)
                    pc.set_clim(-cmax_dZ,cmax_dZ)
                    colorbar(pc, shrink=0.7, label='vertical displacement dz (m)')
                    ax.set_title(f'Vertical deformation dZ at {t} seconds' \
                        + f'\n{event_kinokada}\nKinOkada, without subevents')
                ax.set_xlim(-128.5,-122)
                ax.set_ylim(40,50)
                ax.set_aspect(1/cos(45*pi/180))
                ax.grid(True)
            figs.append(fig)
            close(fig)


        anim = animation_tools.animate_figs(figs, figsize=(10,8))

        # make mp4 file:
        fname = f'{plot_dir}/{event_kinokada}_with_dZ.mp4'
        animation_tools.make_mp4(anim, file_name=fname)
        print('Created ',fname)


    if 1:
        # make animation of slip/vZ (left) with tsunami eta (right)

        cmap_eta = geoplot.tsunami_colormap

        fgno = 1  # which fgout grid
        format = 'binary'  # format of fgout grid output
        outdir = f'{CHT}/geoclaw_runs/offshore_gauges/BR13D_kinokada/_output'
        fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format)
        fgout_grid.read_fgout_grids_data()

        figs = []
        for k,t in enumerate(dtopo.times):
            kv = 2*k  # frame for vZ
            assert waveform_times[kv] == t, \
                    f'*** waveform_time = {waveform_times[kv]} != {t} = t'

            keta = 2*k+1  # fgout frame for tsunami eta plot
            #keta = min(keta, 481)  # fgout ended at t=2400 not 240
            fgout = fgout_grid.read_frame(keta)
            assert fgout.t == t, \
                    f'*** fgout.t = {fgout.t} != {t} = t'

            fig,axs = plt.subplots(1,2,figsize=(10,8))
            for kax,ax in enumerate(axs):
                ax.plot(x_coast,y_coast,'g',linewidth=0.7)
                if kax==0:
                    #pc = ax.pcolormesh(X,Y,seismic_dZ[k,:,:],cmap=cmap)
                    #pc.set_clim(-cmax_dZ,cmax_dZ)
                    #ax.set_title(f'dZ at {t} seconds\n{event}\nSeismic, with subevents')
                    fault.plot_subfaults(axes=ax,slip_color=True,
                                         cmap_slip=cmap_slip,
                                         slip_time=t,plot_box=False,cmin_slip=0.01,
                                         colorbar_shrink=0.7)
                    pc = ax.pcolormesh(X,Y,seismic_vZ[kv,:,:],cmap=cmap_v)
                    pc.set_clim(0,clim_v)
                    ax.set_title(f'Seismic Waves at {t} seconds\n{event}\nWith KinOkada slip')
                if kax==1:
                    # tsunami eta
                    if 0:
                        # plot topo colored by elevation:
                        plottools.pcolorcells(fgout.X,fgout.Y,fgout.B,
                                              ax=ax,cmap=geoplot.land1_colormap)
                        ax.set_clim(0,2000)

                    eta = ma.masked_where(fgout.h<0.001, fgout.eta)

                    eta_plot = plottools.pcolorcells(fgout.X,fgout.Y,eta,
                                                     ax=ax,cmap=cmap_eta)
                    eta_plot.set_clim(-4,4)
                    cb = colorbar(eta_plot, extend='both', shrink=0.7)
                    cb.set_label('Ocean surface elevation (m)')
                    ax.set_title(f'Tsunami Waves at {t} seconds' \
                        + f'\n{event_kinokada}\nKinOkada, without subevents')

                ax.axis('scaled')
                ax.set_xlim(-128.5,-122)
                ax.set_ylim(40,50)
                ax.set_aspect(1/cos(45*pi/180))
                ax.grid(True)
            figs.append(fig)
            close(fig)


        anim = animation_tools.animate_figs(figs, figsize=(10,8))

        # make mp4 file:
        fname = f'{plot_dir}/{event_kinokada}_with_tsunami.mp4'
        animation_tools.make_mp4(anim, file_name=fname)
        print('Created ',fname)


    if 0:
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
            fname = '%s/%s_y%s.png' % (plot_dir,event, str(y0).replace('.','-'))
            savefig(fname)
            print('Created ',fname)
