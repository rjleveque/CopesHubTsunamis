from pylab import *
from clawpack.geoclaw import dtopotools, topotools

dtopo = dtopotools.DTopography('dtopofiles/buried-locking-mur13-deep_instant.dtt3',3)

fname = 'vertical_displacements_sensitivity/vert_displacements_all_xgrid_buried-random-mur13-deep-nosub.dtt3'
dtopo_nosub = dtopotools.DTopography(fname,3)

fname = 'vertical_displacements_sensitivity/vert_displacements_all_xgrid_ft-random-mur13-deep.dtt3'
dtopo_ft = dtopotools.DTopography(fname,3)

def plot_dZ_transect(y0):
    j = where(dtopo.y<y0)[0].max()
    clf()
    plot(dtopo.x,dtopo.dZ[-1,j,:],'b',label='original')
    plot(dtopo.x,dtopo_nosub.dZ[-1,j,:],'r',label='nosub')
    plot(dtopo.x,dtopo_ft.dZ[-1,j,:],'g',label='ft')
    legend()
    grid(True)
    title('dZ transect at latitude y0 = %.2f' % y0);

plot_dZ_transect(45)
savefig('dZy45.png')

plot_dZ_transect(46)
savefig('dZy46.png')
