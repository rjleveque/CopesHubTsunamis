
from pylab import *
import clawpack.pyclaw.gauges as gauges
import os,sys

sys.path.insert(0,'.')
from params import event, location

outdir = os.path.abspath('./_output')
plotdir = os.path.abspath('./_plots')
os.system('mkdir -p %s' % plotdir)
print('Will take output from \n    %sand send plots to \n    %s' \
        % (outdir,plotdir))
        
def make_plot(gaugeno, location, event, outdir, plotdir):
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    x,y = gauge.location
    t = gauge.t / 60.   # convert to minutes
    q = gauge.q
    h = q[0,:]
    #h = where(h>0.01, h, 1.e6)
    u = divide(q[1,:], h, where=h>0.01, out=zeros(h.shape))
    v = divide(q[2,:], h, where=h>0.01, out=zeros(h.shape))
    s = sqrt(u**2 + v**2)
    mflux = h*s
    eta = q[3,:]
    
    maxlevel = gauge.level.max()
    if gauge.level.min() < maxlevel:
        for a in [h,u,v,s,mflux,eta]:
            a = where(gauge.level == maxlevel, a, nan)

    figure(400, figsize=(8,8))
    clf()

    subplot(311)
    plot(t, h, 'b')
    xlabel('')
    ylabel('Flow depth (m)')
    grid(linewidth=0.5)
    title('Gauge %i at x = %.5f, y = %.5f\nIn %s for Event %s' \
            % (gaugeno,x,y,location,event))

    subplot(312)
    plot(t, s, 'b')
    xlabel('')
    ylabel('speed (m/s)')
    grid(linewidth=0.5)
    title('')

    subplot(313)
    plot(t, mflux, 'b')
    xlabel('time (Minutes after earthquake)')
    ylabel('momentum flux (m^3 / s^2)')
    grid(linewidth=0.5)
    title('')

    fname = plotdir + '/%s_%s_Gauge%s.png' \
            % (location,event,str(gaugeno).zfill(5))
    savefig(fname, bbox_inches='tight')
    print('Created %s' % fname)

if __name__ == '__main__':

    gaugenos = range(1001,1051,1)
    
    for gaugeno in gaugenos:
        make_plot(gaugeno, location, event, outdir, plotdir)
