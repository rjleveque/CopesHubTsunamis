"""
Code to make plots for  multiple GeoClaw runs, with a different dtopo file
for each.

It is assumed that geoclaw has already been run using run

It is assumed all parameters set in setrun are the same for all realizations, 
except for the dtopo file to use.

The function make_all_cases_dtopos returns *caselist*, a list of 
dictionaries.  Each dictionary should define whatever parameters
are needed for one case.

The function run_one_case_dtopo(case) takes a dictionary *case* as input
and does what's necessary to run GeoClaw for that one case.

Before running, make sure the following are properly set:

    dtopo_dir is set to the directory containing the dtopo files
    
    runs_dir is set to the destination for _output and _plots directories
        (possibly on a scratch disk).
        
    dtopo_files is set to a list of dtopo files to use
    nprocs is set to the number of geoclaw runs that can be run in parallel.
        (each run will also use OpenMP with OMP_NUM_THREADS threads)
        
Set dry_run = False before executing to actually run GeoClaw.

"""

from numpy import *
import os,sys,glob

dry_run = True  # If True, only print out settings, do not run GeoClaw

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'offshore' # need to change module names below for other locations

sys.path.insert(0, os.path.abspath('..'))
import plot_gauge_max_offshore as plot_gauges
#import process_fgmax_offshore as process_fgmax
import make_fgout_animation_offshore as make_fgout_animation


# location for big files:
this_dir = os.getcwd()
# Randy's laptop:
scratch_dir = this_dir.replace('git/CopesHubTsunamis/geoclaw_runs', \
                               'scratch/CHT_runs')
# for hyak:
scratch_dir = scratch_dir.replace('/mmfs1/home', '/gscratch/tsunami')

# where to find output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)
#runs_dir = os.path.abspath('.')

print('runs_dir = ',runs_dir)

sizes = ['SM','M','L','XL','XXL']
all_events = ['CSZ_%s1_noext' % size for size in sizes] \
           + ['CSZ_%s2_noext' % size for size in sizes] \
           + ['CSZ_%s3_noext' % size for size in sizes]

if 1:
    events = []
    for e in all_events:
        if e[4] in ['M','L']: events.append(e)

if 0:
    # or specify particular events:
    events = ['CSZ_L1_noext']

geoclaw_outputs = os.path.abspath('%s/geoclaw_outputs' % runs_dir)
outdirs = ['%s/_output_%s' % (geoclaw_outputs, event) for event in events]
#print('outdirs = ', outdirs)

geoclaw_plots = os.path.abspath('%s/geoclaw_plots' % runs_dir)
plotdirs = ['%s/_plots_%s' % (geoclaw_plots, event) for event in events]
#plotdirs = [outdir.replace('output','plot') for outdir in outdirs]
print('plotdirs = ', plotdirs)

if dry_run:
    print('DRY RUN - location = %s,  events to process:\n' % location ,events)

def make_html_index(plotdir,event):
    html_fname = os.path.join(plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write('\n<ul>\n<li><a href="fgmax">fgmax plots</a>\n')
        f.write('<li><a href="gauges">gauge plots</a>\n</ul>\n')
    print('Created ',html_fname)


if not dry_run:
    for k in range(len(outdirs)):
        outdir = outdirs[k]
        plotdir = plotdirs[k]
        event = events[k]
        run_name = '%s_%s' % (location,event)

        try:
            if 1:
                gauges_plotdir = plotdir + '/gauges'
                #gaugenos = list(range(1,642,20)) + list(range(2000,2015))
                gaugenos = 'all'
                plot_gauges.make_gauge_plot(gaugenos, outdir, gauges_plotdir,
                                            location, event)

            if 0:
                fgmax_plotdir = plotdir + '/fgmax'
                fg, t_hours = process_fgmax.load_fgmax(outdir)
                process_fgmax.make_fgmax_plots(fg, fgmax_plotdir, run_name, t_hours)

            if 1:
                make_fgout_animation.make_anim(outdir, plotdir, location, event)

            if 0:
                make_html_index(plotdir,event)

        except:
            print('*** Plotting failed for %s'  % run_name)
