
from pylab import *
import clawpack.pyclaw.gauges as gauges
from clawpack.visclaw import gaugetools

def make_gauge_B0(outdir):

    if 1:
        setgauges = gaugetools.read_setgauges(outdir)
        gaugenos = setgauges.gauge_numbers
        print('+++ gaugenos: ', gaugenos)
    else:
        gaugenos = [1001,1002]  # for testing

    fname = 'gauges_B0.csv'
    with open(fname,'w') as f:
        f.write('# gaugeno, longitude, latitude, B0\n')
        for gaugeno in gaugenos:
            gauge = gauges.GaugeSolution(gaugeno, outdir)
            x,y = gauge.location
            t0 = gauge.t[0]
            q = gauge.q
            h = q[0,0]
            eta = q[3,0]
            B0 = eta - h
            level = gauge.level[0]
            print('gaugeno %s at t0 =%8.3f on level %2i has h =%9.3f, eta =%9.3f, B0 =%9.3f' \
                % (gaugeno, t0, level, h, eta, B0))
        
            f.write('%8i,  %14.8f, %14.8f, %10.3f\n' % (gaugeno,x,y,B0))

    print('Created ',fname)


def read_gauge_B0(fname='gauges_B0.csv'):
    
    lines = open(fname).readlines()
    gauge_B0 = {}
    gauge_location = {}
    for line in lines[1:]:
        tokens = line.split(',')
        gaugeno = int(tokens[0])
        x = float(tokens[1])
        y = float(tokens[2])
        B0 = float(tokens[3])
        gauge_B0[gaugeno] = B0
        gauge_location[gaugeno] = (x,y)
    return gauge_B0

