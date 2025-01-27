from pylab import *
import os,sys,glob
from clawpack.geoclaw import dtopotools

#event = 'buried-random-str10-deep'

#depths = ['deep','middle','shallow']
#ruptures = ['buried-locking-str10']

tplot = 180 # seconds after initial rupture time

dtopofiles = glob.glob('dtopofiles_kinematic/*.dtt3')
events = []
for dtopofile in dtopofiles:
    #print('dtopofile: ',dtopofile)
    event = os.path.split(dtopofile)[-1]
    #print('event: ', event)
    event = event.replace('_okada_kinematic.dtt3','')
    #print('event: ', event)
    events.append(event)
events = events[:1]

plotdir = 'plots_compare_seismic_okada_kinematic'
os.system('mkdir -p %s' % plotdir)

for event in events:

    fname = '../dtopofiles/%s.dtt3' % event
    dtopo_seismic = dtopotools.DTopography(fname, 3)

    fname = 'dtopofiles_kinematic/%s_okada_kinematic.dtt3' % event
    dtopo_okk = dtopotools.DTopography(fname, 3)


    fig,(ax0,ax1) = plt.subplots(ncols=2,nrows=1,figsize=(12,9))

    X = dtopo_seismic.X; Y = dtopo_seismic.Y; dZ_at_t = dtopo_seismic.dZ_at_t
    #dz_max = dtopo_seismic.dZ.max()
    #cmax_dZ = dz_max
    cmax_dZ = 10

    #tfinal = dtopo_seismic.times[-1] + 1  # 1 second after final dZ
    dZ_tplot_seismic = dZ_at_t(tplot)
    dz_max = dZ_tplot_seismic.max()
    print('dz_max for seismic at t = %is is %.2f' \
        % (tplot,dz_max))
    dtopotools.plot_dZ_colors(X,Y,dZ_tplot_seismic,axes=ax0,
                              cmax_dZ=cmax_dZ,
                              dZ_interval=200, add_colorbar=True);
    title_text = 'Seafloor deformation at t = %is (seismic)\n%s\nmax dz = %.2fm' \
                    % (tplot,event,dz_max)
    ax0.set_title(title_text)
    ax0.set_ylim(40,50)
    ax0.set_xlim(-128.5,-122)

    X = dtopo_okk.X
    Y = dtopo_okk.Y
    dZ_at_t = dtopo_okk.dZ_at_t
    dz_max = dtopo_okk.dZ.max()
    #tfinal = dtopo_okk.times[-1] + 1  # 1 second after final dZ
    dZ_tplot_okk = dZ_at_t(tplot)
    dz_max = dZ_tplot_okk.max()
    print('dz_max for Okada kinematic at t = %is is %.2f' \
        % (tplot,dz_max))
    dtopotools.plot_dZ_colors(X,Y,dZ_tplot_okk,axes=ax1,
                              cmax_dZ=cmax_dZ,
                              dZ_interval=200, add_colorbar=True);
    title_text = 'Seafloor deformation at t = %is (kinematic Okada)\n%s\nmax dz = %.2fm' \
                    % (tplot,event,dz_max)
    ax1.set_title(title_text)
    ax1.set_ylim(40,50)
    ax1.set_xlim(-128.5,-122)

    fname = '%s/compare_okada_seismic_%s_t%s.png' \
            % (plotdir,event,str(tplot).zfill(3))
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)

    # Transects:

    fig,ax = plt.subplots(ncols=1,nrows=6,sharex=True,figsize=(9,9))

    for k,y0 in enumerate([48.,47,46.,45,44.,43]):
        axk = ax[k]

        dZ = dZ_tplot_seismic
        j = where(dtopo_seismic.y < y0)[0].max()
        axk.plot(dtopo_seismic.x, dZ[j,:], 'b', label='Seismic at y = %.2f' % y0)

        dZ = dZ_tplot_okk
        j = where(dtopo_okk.y < y0)[0].max()
        axk.plot(dtopo_okk.x, dZ[j,:], 'r', label='Okada at y = %.2f' % y0)

        axk.grid(True)
        axk.set_ylim(-5,15)
        axk.set_yticks([-5,0,5,10])
        axk.legend()

    ax[0].set_title('Vertical displacement at t = %is on transects\n' % tplot \
            + event)
    ax[-1].set_xlabel('Longitude x')
    ax[0].set_xlim(-128,-123)
    tight_layout()

    fname = '%s/transects_okada_seismic_kinematic_%s_%s.png' \
            % (plotdir,event,str(tplot).zfill(3))
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)
