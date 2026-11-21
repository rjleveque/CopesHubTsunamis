
from pylab import *
import os, glob
import clawpack.pyclaw.gauges as gauges

outdir = '_output_BL10S'
outdir = '_output'

plotdir = 'gaugeplots_BL10S_seismic'
savefigs = True
if savefigs:
    os.system('mkdir -p %s'  % plotdir)
    print('Plots will be saved to %s' % plotdir)

#gaugenos = range(2001,2010)

gfiles = glob.glob('%s/gauge*.txt' % outdir)
gaugenos = [int(gfile[-9:-4]) for gfile in gfiles]
print('gaugenos = ',gaugenos)

for gaugeno in gaugenos:
    gauge = gauges.GaugeSolution(gaugeno, outdir)

    t = gauge.t / 3600.   # convert to hours
    q = gauge.q
    h = q[0,:]
    eta = q[-1,:]
    fig,ax = subplots(figsize=(12,6))

    offset = h[-1] - eta[-1]
    eta2h = lambda eta: eta + offset
    h2eta = lambda h: h - offset
    ax2 = ax.secondary_yaxis('right', functions=(eta2h,h2eta))

    ax.plot([0,250],[h[0]-offset,h[0]-offset], 'k--', \
            label='initial h = %.2fm' % h[0])

    ax.plot(t, h-offset, 'r', label='water depth h')

    ax.plot(t, eta, 'b', label='surface eta')
    ax.set_ylabel('surface elevation relative to sea level (meters)',fontsize=12)
    ax2.set_ylabel('water depth (meters)',fontsize=12)
    ax.set_xlabel('time (hours)',fontsize=12)

    gauge_string = 'Gauge %s (pnt_%s_%s)' \
                    % (gaugeno, str(gaugeno)[:2],str(gaugeno)[-2:])
    gauge_string = gauge_string + '\nat longitude x = %.5f, latitude y = %.5f' \
                    % (gauge.location[0], gauge.location[1])
    ax.set_title(gauge_string, fontsize=15)

    ax.set_xlim(-0.1, 4.1)

    if 0:
        ylimits = ax.get_ylim()
        ytext = h[0]-offset + 0.005*(ylimits[1]-ylimits[0])
        ax.text(3.1, ytext, 'initial water depth = %.2fm' % h[0])


    ax.legend(loc='upper right', framealpha=1)

    ax.grid(True)

    if savefigs:
        fname = '%s/Gauge%s.png' % (plotdir,str(gaugeno).zfill(4))
        savefig(fname, bbox_inches='tight')
        print('Created %s' % fname)
        close(fig)

