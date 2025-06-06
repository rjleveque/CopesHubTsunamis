"""
setplot function with case argument for looping over many events from
run_geoclaw_make_plots_dtopos.py

This version makes no frame plots, only fgmax and gauge plots!
"""

import os,sys

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'offshore_gauges' # need to change module names below for other locations

sys.path.insert(0, os.path.abspath('..'))
#import plot_gauge_max_offshore as plot_gauges
import plot_fgmax

# for FrontalThrust (modify max amplitude?):
#import plot_gauge_max_offshore_ft as plot_gauges
import plot_gauges

#import process_fgmax_offshore as process_fgmax
import make_fgout_animation_offshore as make_fgout_animation


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

    try:
        if 1:
            gauges_plotdir = plotdir + '/gauges'
            #gaugenos = list(range(1,642,20)) + list(range(2000,2015))
            gaugenos = 'all'
            plot_gauges.make_gauge_plot(gaugenos, outdir, gauges_plotdir,
                                        location, event)

        if 1:
            fgmax_plotdir = plotdir + '/fgmax'
            #fg, t_hours = process_fgmax.load_fgmax(outdir)
            #process_fgmax.make_fgmax_plots(fg, fgmax_plotdir, run_name, t_hours)
            plot_fgmax.make_fgmax_plot(outdir, fgmax_plotdir,
                                       location, event, dtopofile)

        if 1:
            try:
                make_fgout_animation.make_anim(outdir, plotdir, location, event)
            except:
                print('*** filed to make fgout animation')

        if 1:
            make_html_index(plotdir,event)

    except:
        print('*** Plotting failed for %s'  % run_name)
        raise()

    plotdata = None

    # add code to set plotdata (as in other setplot.py modules) if you also
    # want to make time frame plots.

    return plotdata

def make_html_index(plotdir,event):
    html_fname = os.path.join(plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write('\n<ul>\n<li><a href="fgmax">fgmax plots</a>\n')
        f.write('<li><a href="gauges">gauge plots</a>\n</ul>\n')
    print('Created ',html_fname)

