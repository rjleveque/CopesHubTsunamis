

import sys, os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *

import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import topotools
from clawpack.visclaw import plottools, geoplot, gaugetools
from clawpack.geoclaw import fgout_tools

#event = os.path.split(os.getcwd())[-1]
#event = 'Comparison of different events'
location = 'Seaside'

all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


if 0:

    models = all_models[0]
    subtitle = 'buried-random models'
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

if 1:

    models = all_models[:1]
    subtitle = 'buried-random models'
    events = ['%s-deep' % model for model in models]

print('Events: ', events)

#color_list = 3*['r','b','g','r','b','g','r','b','g']
#linestyle_list = 3*['-','-','-','--','--','--','-.','-.','-.']

plotdir = '.'

def read_gauge(outdir, gaugeno):
    from numpy import divide
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    t = gauge.t
    eta = gauge.q[3,:]
    h = gauge.q[0,:]
    hu = gauge.q[1,:]
    hv = gauge.q[2,:]
    u = divide(hu, h, where=h>0, out=zeros(hu.shape))
    v = divide(hv, h, where=h>0, out=zeros(hv.shape))
    s = sqrt(u**2 + v**2)
    xg = gauge.location[0]
    yg = gauge.location[1]
    return t, eta, s, xg, yg
    
gaugenos = range(1,5,1)

for event in events:
    print('Event ',event)
    event_instant = event + '_instant'
    outdir_event = '/gscratch/tsunami/rjl/CopesHubTsunamis/geoclaw_runs/sites/seaside/multirun2/geoclaw_plots/_plots_%s' \
        % event
    outdir_instant = outdir_event + '_instant'
    
    # test:
    outdir_event = '/Users/rjl/scratch/CHT_runs/sites/seaside/multirun/geoclaw_outputs/_output_buried-random-str10-middle'
    outdir_instant = '/Users/rjl/scratch/CHT_runs/sites/seaside/multirun/geoclaw_outputs/_output_buried-random-str10-shallow'
    gaugenos = [1001]

    for gaugeno in gaugenos:
        figure(205,figsize=(12,8))
        clf()

        t,eta,s,xg,yg = read_gauge(outdir_event, gaugeno)

        subplot(211)
        plot(t/60., eta, color='r', linestyle='-', label='kinematic')
        
        subplot(212)
        plot(t/60., s, color='r', linestyle='-', label='kinematic')
        
        t_instant,eta_instant,s_instant,xg,yg = \
                                read_gauge(outdir_instant, gaugeno)
        
        eta_instant[0] = eta[0]  # preseismic
        
        subplot(211)
        plot(t_instant/60., eta_instant, color='r', linestyle='--', label='instant')
        grid(True)
        xlim(-1,60)
        #xlabel('time (minutes)')
        ylabel('surface elevation (m)')
        legend(framealpha=1)
        title('%s Gauge %s at (%.5f, %.5f)\n%s' \
                % (location,gaugeno,xg,yg,event))

        subplot(212)
        plot(t_instant/60., s_instant, color='r', linestyle='--', label='instant')
        grid(True)
        xlim(-1,60)
        xlabel('time (minutes)')
        ylabel('speed (m/s)')
        legend(framealpha=1)
        
        tight_layout()
        
        fname = '%s/gauge%s_comparison.png' \
                    % (plotdir,str(gaugeno).zfill(5))
        savefig(fname)
        print('Created ',fname)
