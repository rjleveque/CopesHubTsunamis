"""
Sample post-processing script
In practice replace this with a script to plot gauges, fgmax, etc.
"""

from pylab import *
import clawpack.pyclaw.gauges as gauges


def test(outdir, plotdir, location, event):
    """
    Write a simple text file to check that it shows up in plotdir as expected
    Create a simple gauge plot
    """

    outfile = f'{plotdir}/post-process-info.txt'
    with open(outfile, 'w') as f:
        f.write('Sample post-processing output for:\n')
        f.write(f'  event = {event} at location = {location}\n')
        f.write(f'output can be found in outdir:\n    {outdir}\n')
        f.write(f'plots can be found in plotdir:\n    {plotdir}\n')

    print(f'Sample post-processing file written to {outfile}')

    # plot gauge output:
    gaugeno = 101
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    eta = gauge.q[-1,:]
    plot(gauge.t, eta, 'b')
    xlabel('time (seconds)')
    ylabel('eta (m)')
    title(f'Gauge {gaugeno} for event {event}')
    
    fname = f'{plotdir}/{event}_gauge{gaugeno}.png'
    savefig(fname)
    print(f'Created gauge plot {fname}')
