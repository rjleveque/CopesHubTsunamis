"""
Used for AGU24 poster results
"""

from pylab import *
import clawpack.pyclaw.gauges as gauges
import os,sys

        
def make_txt(gaugeno, location, event, outdir, plotdir):
    
    os.system('mkdir -p %s' % plotdir)
    print('Will take output from \n    %sand send plots to \n    %s' \
            % (outdir,plotdir))
            
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    x,y = gauge.location
    #t = gauge.t / 60.   # convert to minutes
    t = gauge.t  # seconds
    q = gauge.q
    h = q[0,:]
    u = divide(q[1,:], h, where=h>0.01, out=zeros(h.shape))
    v = divide(q[2,:], h, where=h>0.01, out=zeros(h.shape))
    #s = sqrt(u**2 + v**2)
    #mflux = h*s
    eta = q[-1,:]
    
    maxlevel = gauge.level.max()
    if gauge.level.min() < maxlevel:
        for a in [h,u,v,s,mflux,eta]:
            a = where(gauge.level == maxlevel, a, nan)

    gauge_out = vstack((t,h,u,v,eta)).T
    fname_gauge_out = '%s/%s_Gauge%s.txt' % (plotdir,event,gaugeno)
    header = 'Event = %s\nGauge %i,  location = %11.6f, %11.6f\n' \
                % (event, gaugeno, x,y) \
             + 't(sec),   h(m),   u(m/s),   v(m/s),   eta(m)'
    savetxt(fname_gauge_out, gauge_out, fmt='%10.3f', header=header)
    
    print('Created %s' % fname_gauge_out)

if __name__ == '__main__':

    sys.path.insert(0,'.')
    from params import event, location

    outdir = os.path.abspath('./_output')
    #outdir = os.path.abspath('./_output_manning025') # for buried-random-mur13-deep
    plotdir = os.path.abspath('./_plots')
    
    gaugenos = range(1001,1051,1)
    gaugenos = [1015, 1045]
    
    for gaugeno in gaugenos:
        make_txt(gaugeno, location, event, outdir, plotdir)

    close('all')
