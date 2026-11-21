from pylab import *
from clawpack.pyclaw.gauges import GaugeSolution

event = 'BL10D'

outdir1 = f'./_output_{event}_seismic'
outdir2 = f'./_output_{event}'
#gaugeno1 = 101

for gaugeno1 in [101,102,103]:
    gaugeno2 = gaugeno1

    gauge1 = GaugeSolution(gauge_id=gaugeno1, path=outdir1)
    gauge2 = GaugeSolution(gauge_id=gaugeno2, path=outdir2)

    figure(500,figsize=(10,7))
    clf()
    subplot(211)

    plot(gauge1.t/3600., gauge1.q[-1,:], 'r', label='Seismic')
    plot(gauge2.t/3600., gauge2.q[-1,:], 'b', label='No subevents, KinOkada')

    #legend(loc='upper right')
    legend()
    xlim(0,2)
    #ylim(-2,4)
    grid(True)
    #xlabel('Hours after earthquake')
    ylabel('meters')

    title(f'Surface eta at Gauge {gaugeno1}')

    subplot(212)

    plot(gauge1.t/3600., gauge1.q[0,:], 'r', label='Seismic')
    plot(gauge2.t/3600., gauge2.q[0,:], 'b', label='No subevents, KinOkada')

    #legend(loc='upper right')
    legend()
    xlim(0,2)
    #ylim(-2,4)
    grid(True)
    xlabel('Hours after earthquake')
    ylabel('meters')

    title(f'Water depth at Gauge {gaugeno1}')

    tight_layout()

    fname = f'Gauge{gaugeno1}_SeismicVsNosub_{event}.png'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)