"""
setplot function with case argument for looping over many events from
run_geoclaw_make_plots_dtopos.py

This version makes no frame plots, only does post-processing
by calling other script(s).
"""

import os,sys
#import plot_gaugereport   # Loyces gauge report
import plot_play_gaugereport   # Loyces gauge report

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'AberHospital' # set to site location if desired to use in plots/filenames


def setplot(plotdata, case):

    """
    Run the code for post-processing.
    This version does not return plotdata for making time frame plots.
    You can also add code to set plotdata and return it if desired.

    case is a dictionary that will be passed in, and it is assumed that
    this dictionary has at least the keys 'outdir', 'plotdir', 'dtopofile'
    with the appropriate values for one particular case (dtopofile).
    """

    from pathlib import Path

    # these parameters should be set properly when calling this setplot:

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofile = case['dtopofile']

    # extract event name:
    if dtopofile is not None:
        event = Path(dtopofile).stem  # drop path and .dtt3 extension
    else:
        event = 'NO_DTOPO'

    run_name = '%s_%s' % (location,event)
    print('In setplot: run_name = ',run_name)

    # call post-processing script(s). This prints a .csv and a .txt file:
    #plot_gaugereport.report(outdir, plotdir, location, event, dtopofile, run_name)
    plot_play_gaugereport.report(outdir, plotdir, location, event, dtopofile, run_name)

    plotdata = None

    # add code to set plotdata (as in other setplot.py modules) if you also
    # want to make time frame plots.

    return plotdata
