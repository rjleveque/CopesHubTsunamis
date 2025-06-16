"""
Code to perform multiple GeoClaw runs, with a different dtopo file
for each, and then make plots for each run.

If run_code is False the code will not be run (e.g. if the runs have been
done and only plotting should be done).

If make_plots is False the plotting will not be done.

It is assumed all parameters set in setrun are the same for all realizations,
except for the dtopo file to use.  The setrun_case.py module contains a
function setrun that takes case as a parameter so that case['dtopofiles']
can be used in setrun.

For making plots, it is assumed a setplot_case.py module contains a
function setplot that takes case as a parameter so that the event name
can be extracted for adding to plots / filenames.

The function make_all_cases_dtopos returns *caselist*, a list of
dictionaries.  Each dictionary should define whatever parameters
are needed for one case.

The function run_one_case_clawpack(case) takes a dictionary *case* as input
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

# what to do:
run_code = True
make_plots = True

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.join(root_dir, 'common_code'))
from cases_dtopos import make_all_cases_dtopos

if 0:
    # if clawmultip_tools and multip_tools should import from clawpack:
    sys.path.insert(0, '/Users/rjl/git/clawpack/clawmultip/src/python/clawmultip')
    # (Otherwise use the versions in common_code)

from clawmultip_tools import run_one_case_clawpack
from multip_tools import run_many_cases_pool

# location for big files:
this_dir = os.getcwd()
# Randy's laptop:
scratch_dir = this_dir.replace('git/CopesHubTsunamis/geoclaw_runs', \
                               'scratch/CHT_runs')
#scratch_dir = '/Users/rjl/tests/CHT_runs'

# for hyak:
scratch_dir = scratch_dir.replace('/mmfs1/home', '/gscratch/tsunami')

# where to find all the dtopo files:
dtopo_dir = os.path.join(root_dir, \
                       'dtopo/CSZ_groundmotions/dtopofiles_3D_5km')

# for hyak:
dtopo_dir = dtopo_dir.replace('/mmfs1/home', '/gscratch/tsunami')
print('dtopo_dir = ',dtopo_dir)

# where to put output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)
#runs_dir = os.path.abspath('.')

# create scratch directory for output, if it doesn't exist:
os.system('mkdir -p %s' % runs_dir)

# path to geoclaw executable:
# (should agree with how EXE is set in Makefile)
if run_code:
    xgeoclaw_path = os.path.join(root_dir, 'geoclaw_runs/xgeoclaw_v512')
else:
    xgeoclaw_path = None  # do not run GeoClaw code

# number of events to run and/or plot simultaneously:
#nprocs = 9  # now set by reading sys.argv in __main__

# Specify the list of events to loop over for geoclaw runs:



if 1:

    all_models = \
        ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
         'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


    models = all_models
    #models = all_models[:3]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

    events.sort()

    #events = events[:9]

    events = [e.replace('-','_') for e in events]
    events = [e+'_3D_5km' for e in events]

    instant = True
    if instant:
        events = [e+'_instant' for e in events]



dtopo_files = ['%s/%s.dtt3' % (dtopo_dir,f) for f in events]

#print('dtopo_files = ',dtopo_files)
#print('events = ',events)

if __name__ == '__main__':

    import sys
    try:
        nprocs = int(sys.argv[1])
    except:
        raise Exception('*** Missing integer argument nprocs on command line')


    print('\n--------------------------')
    if dry_run:
        # just print out settings, no runs...
        print('DRY RUN - settings in run_geoclaw_dtopos.py')

    print('Will run GeoClaw for %i dtopo files' % len(dtopo_files))
    print('dtopo files should be in dtopo_dir:\n    ', dtopo_dir)
    print('list of dtopo_files to process:')
    for fpath in dtopo_files:
        print('  %s' %  os.path.split(fpath)[-1])
        if not os.path.isfile(fpath):
            print('    *** file not found')
    print('output will go in \n    %s/geoclaw_outputs/' % runs_dir)
    print('nprocs = %i jobs will run simultaneously' % nprocs)
    print('GeoClaw executable:\n    ',xgeoclaw_path)

    if dry_run:
        print('Set dry_run=False and re-execute to run GeoClaw')
        print('--------------------------')
    else:
        # make list of dictionaries with parameters for each case:
        caselist = make_all_cases_dtopos(dtopo_dir, dtopo_files, runs_dir,
                                         xgeoclaw_path, make_plots)

        # run all cases using nprocs processors:
        run_many_cases_pool(caselist, nprocs, run_one_case_clawpack)
