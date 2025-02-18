"""
setplot function with case argument for looping over many events from
run_geoclaw_make_plots_dtopos.py
"""

import os,sys

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'Newport' # need to change module names below for other locations

sys.path.insert(0, root_dir + '/common_code')
import plot_gauges_site  # now works for all sites

sys.path.insert(0, os.path.abspath('..'))
import process_fgmax_Newport as process_fgmax
import make_fgout_animation_Newport as make_fgout_animation


def setplot(plotdata, case={}):

    """
    Run the code for making fgmax and gauge plots and fgout animations.
    This version does not return plotdata for making time frame plots.
    Add code to set plotdata and return it if desired.
    """

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofiles = case['dtopofiles']

    dtopofile = dtopofiles[0][-1]
    event = os.path.splitext(os.path.split(dtopofile)[-1])[0]

    run_name = '%s_%s' % (location,event)
    print('In setplot: run_name = ',run_name)

    if 1:
        gauges_plotdir = plotdir + '/gauges'
        os.system('mkdir -p %s' % gauges_plotdir)
        plot_gauges_site.make_all_plots_and_report(outdir, gauges_plotdir,
                         location=location, event=event,
                         gaugenos='all', sea_level=0.)

    if 1:
        try:
            #fgmax_plotdir = plotdir + '/fgmax'
            fgmax_plotdir = plotdir  # fgmax1 or fgmax2 now added in function
            os.system('mkdir -p %s' % fgmax_plotdir)
            process_fgmax.make_all_fgmax_plots(outdir, fgmax_plotdir,
                                           location=location, event=event)
        except:
            print('*** problem with process_fgmax in %s, skipping' \
                  % run_name)
    if 1:
        make_fgout_animation.make_anim(outdir, plotdir, location, event)

    if 0:
        make_html_index(plotdir,event)

    plotdata = None

    # add code to set plotdata (as in other setplot.py modules) if you also
    # want to make time frame plots.

    return plotdata
