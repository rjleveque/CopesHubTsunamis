"""
setplot function with case argument for looping over many events from
run_geoclaw_make_plots_dtopos.py

This version makes no frame plots, only does post-processing
by calling other script(s).
"""

import os,sys
import post_process   # sample script for illustration

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'test' # set to site location if desired to use in plots/filenames


def setplot(plotdata, case={}):

    """
    Run the code for post-processing.
    This version does not return plotdata for making time frame plots.
    You can also add code to set plotdata and return it if desired.
    """

    # these parameters should be set by runclaw_makeplots_dtopos.py
    # when calling this setplot:

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofiles = case['dtopofiles']

    # extract event name:
    dtopofile = dtopofiles[0][-1]
    event = os.path.splitext(os.path.split(dtopofile)[-1])[0]

    run_name = '%s_%s' % (location,event)
    print('In setplot: run_name = ',run_name)

    # call post-processing script(s). This one just prints a text file:
    post_process.test(outdir, plotdir, location, event)

    plotdata = None

    # add code to set plotdata (as in other setplot.py modules) if you also
    # want to make time frame plots.

    return plotdata

