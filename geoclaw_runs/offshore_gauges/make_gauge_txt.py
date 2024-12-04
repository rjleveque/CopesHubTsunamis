from pylab import *
import os

#event = os.path.split(os.getcwd())[-1]
event = 'buried-random-mur13-deep'

# make .txt file with gauge maxima for Yong's AGU poster:

outdir1 = 'buried-random-mur13-deep_shallow/_output_15sec'
outdir2 = 'buried-random-mur13-deep_instant/_output'
gaugenos = range(1,660,1)

plotdir = '_plots_agu24a'
os.system('mkdir -p %s' % plotdir)

def read_gauges(outdir, gaugenos=None):
    from clawpack.pyclaw import gauges
    from clawpack.visclaw import gaugetools
    if gaugenos is None:
        setgauges = gaugetools.read_setgauges(outdir)
        gaugenos = setgauges.gauge_numbers
    
    xg = []
    yg = []
    gmax = []
    for gaugeno in gaugenos:
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        h = gauge.q[0,:]
        eta = where(h>0, gauge.q[-1,:], nan)
        gmax.append(nanmax(eta))
        #gmax.append(eta.max())
        xg.append(gauge.location[0])
        yg.append(gauge.location[1])

    # sort gauges from south to north
    isort = argsort(yg)
    gaugenos = [gaugenos[i] for i in isort]
    #print('+++ gaugenos sorted S to N:\n  ',gaugenos)
    gmax = array(gmax)
    xg = array(xg)
    yg = array(yg)
    xg = xg[isort]
    yg = yg[isort]
    gmax = gmax[isort]

    return xg,yg,gmax,gaugenos
    
gauge_max_eta = None
for k,outdir in enumerate([outdir1,outdir2]): 
    print('Reading gauges from ',outdir)
    xg,yg,gmax,gaugenos = read_gauges(outdir,gaugenos)             
    if gauge_max_eta is None:
        gauge_max_eta = empty((len(yg), 4), dtype=float)
        gauge_max_eta[:,0] = xg
        gauge_max_eta[:,1] = yg
    gauge_max_eta[:,2+k] = gmax
    
fname = '%s/GeoClawOffshoreGauges_%s.txt' % (plotdir,event)
header = 'Event %s non-dispersive GeoClaw with 15 arcsec resolution\n' % event \
            + 'Lon  Lat  3D_kinematic  3D_static'
savetxt(fname, gauge_max_eta, fmt='%13.6f', header=header)
print('Created ',fname)
