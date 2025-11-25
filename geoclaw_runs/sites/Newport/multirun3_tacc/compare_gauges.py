"""
Script to compare time series at gauges between two different GeoClaw runs.
with results in outdir and outdir2.

Modify outdir and outdir2 to compare runs varying different things.
You might also want to modify `labels` to give properly descriptive labels
in the figure legends.
"""

from pylab import *
import os
from clawpack.pyclaw.gauges import GaugeSolution

event = 'BL10D'

outdir = f'geoclaw_outputs/_output_{event}'

# compare to:
#outdir2 = '../multirun1/geoclaw_outputs/_gauges_buried-locking-mur13-deep'
#outdir2 = '../multirun1/geoclaw_outputs/_output_buried-locking-str10-deep'
outdir2 = '../multirun2_tacc/geoclaw_outputs/_output_{event}'


outdirs = [outdir, outdir2]
c = ['b','r']  # colors for lines
labels = ['TACC', 'Hyak']


gaugenos = [1007,1022,1041]

for kg,gaugeno in enumerate(gaugenos):
    figure(gaugeno, figsize=(10,8))
    clf()
    for ke,outdir in enumerate(outdirs):
        gauge = GaugeSolution(gauge_id=gaugeno, path=outdir)
        exlabel = labels[ke]

        maxlevel = gauge.level.max()
        print(f'Filtering out where gauge.level < {maxlevel}')
        i = where(gauge.level < maxlevel)
        gauge.q[:,i] = nan

        subplot(211)
        plot(gauge.t/60., gauge.q[-1,:], c[ke], label=exlabel)
        title_string = f'Surface eta at Gauge {gaugeno}'
        title(title_string, fontsize=15)
        grid(True)
        legend(loc='upper right', framealpha=1)
        ylabel('meters')

        subplot(212)
        plot(gauge.t/60., gauge.q[0,:], c[ke], label=exlabel)
        title_string = f'Water depth at Gauge {gaugeno}'
        title(title_string, fontsize=15)
        grid(True)
        legend(loc='upper right', framealpha=1)
        ylabel('meters')

    xlabel('Minutes after earthquake')

    tight_layout()

    if 1:
        fname = f'Gauge{gaugeno}_Comparison_{event}.png'
        savefig(fname)
        print('Created ', fname)
