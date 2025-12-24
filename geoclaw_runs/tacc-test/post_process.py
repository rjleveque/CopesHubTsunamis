"""
Sample post-processing script
In practice replace this with a script to plot gauges, fgmax, etc.
"""

def test(outdir, plotdir, location, event):
    """
    write a simple text file to check that it shows up in plotdir as expected
    """

    outfile = f'{plotdir}/post-process-info.txt'
    with open(outfile, 'w') as f:
        f.write('Sample post-processing output for:\n')
        f.write(f'  event = {event} at location = {location}\n')
        f.write(f'output can be found in outdir:\n    {outdir}\n')
        f.write(f'plots can be found in outdir:\n    {plotdir}\n')

    print(f'Sample post-processing file written to {outfile}')
