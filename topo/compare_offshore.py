"""
Compare several topo DEMs
plot_topo_transect function plots comparison on same plot
"""

from pylab import *
from clawpack.geoclaw import topotools, dtopotools, kmltools

if 0:
    gebco = topotools.Topography(\
        '/Users/rjl/git/CopesHubTsunamis/topo/topofiles/gebco_2021_n54.0_s38.0_w-137.0_e-122.0.asc',3)
    gebcoc = gebco.crop([-127,-123,45,49])

    etopo = topotools.Topography(\
        '/Users/rjl/topo/topofiles/etopo1_-163_-122_38_63.asc',3)
    etopoc = etopo.crop([-127,-123,45,49])

    topo3s = topotools.Topography(\
        '/Users/rjl/topo/topofiles/nw_pacific_3sec_SW_WA.asc',3)

    #etopo22 = topotools.Topography('topofiles/etopo22a.asc',3)

    etopo22 = topotools.Topography(\
        'topofiles/etopo22_-130_-123_38_52.asc',3)

    etopo22_1min = topotools.Topography(\
        'topofiles/etopo22_1min_-130_-123_38_52.asc',3)

y0 = 46.25

def plot_topo_transect(y0):
    je = where(etopoc.y<y0)[0].max()
    jg = where(gebcoc.y<y0)[0].max()
    j3 = where(topo3s.y<y0)[0].max()
    j22 = where(etopo22.y<y0)[0].max()
    j221m = where(etopo22_1min.y<y0)[0].max()
    
    figure(figsize=(10,5))
    plot(etopoc.x,etopoc.Z[je,:],'b',label='etopo1 60 arcsec')
    plot(etopo22_1min.x,etopo22_1min.Z[j221m,:],'r',label='etopo22 60 arcsec')
    plot(etopo22.x,etopo22.Z[j22,:],'k',label='etopo22 15 arcsec')
    plot(gebcoc.x,gebcoc.Z[jg,:],'g',label='gebco 15 arcsec')
    #plot(topo3s.x,topo3s.Z[j3,:],'m',label='coastal 3 arcsec')
    
    grid(True)
    xlim(-126.5, -123.5)
    legend(loc='upper left')
    title('Topo transects at y = %.4f' % y0, fontsize=15)
    xlabel('longitude', fontsize=12)
    ylabel('meters', fontsize=12)
