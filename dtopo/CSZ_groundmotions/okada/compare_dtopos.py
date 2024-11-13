from pylab import *
from clawpack.geoclaw import dtopotools

event = 'buried-random-str10-deep'

fname = '../dtopofiles/%s_instant.dtt3' % event
dtopo_instant = dtopotools.DTopography(fname, 3)

fname = '%s_okada_instant.dtt3' % event
dtopo_oki = dtopotools.DTopography(fname, 3)


fig,(ax0,ax1) = plt.subplots(ncols=2,nrows=1,figsize=(12,9))

X = dtopo_instant.X; Y = dtopo_instant.Y; dZ_at_t = dtopo_instant.dZ_at_t
dz_max = dtopo_instant.dZ.max()
#cmax_dZ = dz_max
cmax_dZ = 15
print('dz_max for gm instant is %.2f, using cmax_dZ = %.2f' % (dz_max,cmax_dZ))

tfinal = dtopo_instant.times[-1] + 1  # 1 second after final dZ
dtopotools.plot_dZ_colors(X,Y,dZ_at_t(tfinal),axes=ax0, 
                          cmax_dZ=cmax_dZ, 
                          dZ_interval=200, add_colorbar=True);
title_text = 'Seafloor deformation (static seismic)\n%s\nmax dz = %.2fm' \
                % (event,dz_max)
ax0.set_title(title_text)
ax0.set_ylim(40,50)
ax0.set_xlim(-128.5,-122)

X = dtopo_oki.X
Y = dtopo_oki.Y
dZ_at_t = dtopo_oki.dZ_at_t
dz_max = dtopo_oki.dZ.max()
tfinal = dtopo_oki.times[-1] + 1  # 1 second after final dZ
dtopotools.plot_dZ_colors(X,Y,dZ_at_t(tfinal),axes=ax1, 
                          cmax_dZ=cmax_dZ, 
                          dZ_interval=200, add_colorbar=True);
title_text = 'Seafloor deformation (static Okada)\n%s\nmax dz = %.2fm' \
                % (event,dz_max)
ax1.set_title(title_text)
ax1.set_ylim(40,50)
ax1.set_xlim(-128.5,-122)

fname = 'compare_okada_seismic_instant_%s.png' % event
savefig(fname, bbox_inches='tight')
print('Created ', fname)

# Transects:

fig,ax = plt.subplots(ncols=1,nrows=6,sharex=True,figsize=(9,9))

for k,y0 in enumerate([48.,47,46.,45,44.,43]):
    axk = ax[k]
    
    dtopo = dtopo_instant
    j = where(dtopo.y < y0)[0].max()
    axk.plot(dtopo.x, dtopo.dZ[-1,j,:], 'b', label='Seismic at y = %.2f' % y0)
    
    dtopo = dtopo_oki
    j = where(dtopo.y < y0)[0].max()
    axk.plot(dtopo.x, dtopo.dZ[-1,j,:], 'r', label='Okada at y = %.2f' % y0)
    
    axk.grid(True)
    axk.set_ylim(-5,15)
    axk.set_yticks([-5,0,5,10])
    axk.legend()
    
ax[0].set_title('Vertical (static) displacement on transects\n' +\
        event)
ax[-1].set_xlabel('Longitude x')
ax[0].set_xlim(-128,-123)
tight_layout()

fname = 'transects_okada_seismic_instant_%s.png' % event
savefig(fname, bbox_inches='tight')
print('Created ', fname)
