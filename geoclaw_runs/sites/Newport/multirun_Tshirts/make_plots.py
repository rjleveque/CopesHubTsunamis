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

#dry_run = True  # If True, only print out settings, do not run GeoClaw
dry_run = False  # If True, only print out settings, do not run GeoClaw

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

location = 'Newport' # need to change module names below for other locations

sys.path.insert(0, root_dir + '/common_code')
import plot_gauges_site  # now works for all sites

sys.path.insert(0, os.path.abspath('..'))
import process_fgmax_Newport as process_fgmax
import make_fgout_animation_Newport as make_fgout_animation


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

# where to find all the dtopo files:
dtopo_dir = os.path.join(root_dir, 'dtopo/CSZ_Bandon_noext/')

# for hyak:
dtopo_dir = dtopo_dir.replace('/mmfs1/home', '/gscratch/tsunami')
print('dtopo_dir = ',dtopo_dir)

sizes = ['SM','M','L','XL','XXL']

all_events = []
for size in sizes:
    for M in [1,2,3]:
        all_events.append('CSZ_%s%s_noext' % (size,M))

if 0:
    events = all_events

if 1:
    events = []
    for e in all_events:
        if e[4] in ['M','L']: events.append(e)

if 0:
    # or specify particular events:
    events = ['CSZ_L1_noext','CSZ_XL1_noext']

dtopo_files = ['%s/%s.tt3' % (dtopo_dir,f) for f in events]


geoclaw_outputs = os.path.abspath('%s/geoclaw_outputs' % runs_dir)
outdirs = ['%s/_output_%s' % (geoclaw_outputs, event) for event in events]
#print('outdirs = ', outdirs)

geoclaw_plots = os.path.abspath('%s/geoclaw_plots' % runs_dir)
plotdirs = ['%s/_plots_%s' % (geoclaw_plots, event) for event in events]
#plotdirs = [outdir.replace('output','plot') for outdir in outdirs]
print('plotdirs = ', plotdirs)

#gaugenos = range(1001,1079,1)
#print('Will make %i gauge plots for each event' % len(gaugenos))

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

        if 1:
            gauges_plotdir = plotdir + '/gauges'
            os.system('mkdir -p %s' % gauges_plotdir)
            plot_gauges_site.make_all_plots_and_report(outdir, gauges_plotdir, 
                             location=location, event=event,
                             gaugenos='all', sea_level=0.)

        if 1:
            try:
                fgmax_plotdir = plotdir + '/fgmax'
                process_fgmax.make_all_fgmax_plots(outdir, fgmax_plotdir,
                                               location=location, event=event)
            except:
                print('*** problem with process_fgmax in %s, skipping' \
                      % run_name)

        if 1:
            make_fgout_animation.make_anim(outdir, plotdir, location, event)

        if 0:
            make_html_index(plotdir,event)
