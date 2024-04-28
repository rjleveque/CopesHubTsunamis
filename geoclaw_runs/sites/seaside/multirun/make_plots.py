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

dry_run = False  # If True, only print out settings, do not run GeoClaw

location = 'Seaside'

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.abspath('..'))
import plot_gauges, process_fgmax

# location for big files:
this_dir = os.getcwd()
# Randy's laptop:
scratch_dir = this_dir.replace('git/CopesHubTsunamis/geoclaw_runs', \
                               'scratch/CHT_runs')
# for hyak:
scratch_dir = this_dir.replace('/mmfs1/home', '/gscratch/tsunami')

# where to find output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)
#runs_dir = os.path.abspath('.')

#events = ['buried-random-str10-middle','buried-random-str10-shallow']

events = ['buried-locking-mur13-deep', 'buried-locking-skl16-deep', 
          'buried-locking-mur13-middle', 'buried-locking-skl16-middle',
          'buried-locking-mur13-shallow']

output_dirs = os.path.abspath('%s/geoclaw_outputs' % runs_dir)
outdirs = ['%s/_output_%s' % (output_dirs, event) for event in events]
#print('outdirs = ', outdirs)

plot_dirs = os.path.abspath('%s/geoclaw_plots' % runs_dir)
plotdirs = ['%s/_plots_%s' % (plot_dirs, event) for event in events]
#plotdirs = [outdir.replace('output','plot') for outdir in outdirs]
#print('plotdirs = ', plotdirs)

gaugenos = range(1001,1051,1)
print('Will make %i gauge plots for each event' % len(gaugenos))

if dry_run:
    print('DRY RUN - location = %s,  events to process:\n' % location ,events)

if not dry_run:
    for k in range(len(outdirs)):
        outdir = outdirs[k]
        plotdir = plotdirs[k]
        event = events[k]
        run_name = '%s_%s' % (location,event)

        gauges_plotdir = plotdir + '/gauges'
        for gaugeno in gaugenos:
            plot_gauges.make_plot(gaugeno, location, event, outdir,
                                  gauges_plotdir)

        fgmax_plotdir = plotdir + '/fgmax'
        fg, t_hours = process_fgmax.load_fgmax(outdir)
        process_fgmax.make_fgmax_plots(fg, fgmax_plotdir, run_name, t_hours)
        process_fgmax.make_kmz_plots(fg, fgmax_plotdir, run_name)
