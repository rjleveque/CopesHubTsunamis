"""
Code to perform multiple GeoClaw runs, with a different dtopo file
for each, and then make plots for each run.

If run_code is False the code will not be run (e.g. if the runs have been
done and only plotting should be done).

If make_plots is False the plotting will not be done.

It is assumed all parameters set in setrun are the same for all realizations,
except for the dtopo file to use.  The setrun_case.py module contains a
function setrun that takes case as a parameter so that case['dtopofile']
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

Note that the call to multip_tools.run_many_cases_pool should be in __main__

This script can be executed from the command line as:
    python runclaw_makeplots_dtopos.py NPROCS
where NPROCS is the number of jobs to run in parallel using the
Python mulitprocessing tools.

Alternatively, 
    python runclaw_makeplots_dtopos.py NPROCS FIRST_EVENT LAST_EVENT
will restrict the set of events / dtopo_files created to this subset.

In this script `all_events` is set to be the set of 36 CoPes Hub ground motions
For the naming and numbering convention, see
    https://depts.washington.edu/ptha/CHTuser/dtopo/ground-motions/
The integers FIRST_EVENT  and LAST_EVENT can be used to specify a subset.
This is useful for starting multiple jobs via slurm to run them all.

For example, 
    python runclaw_makeplots_dtopos.py 2 1 4
will run 2 geoclaw jobs at a time in parallel on the first 4 events
   1    BL10D    buried-locking-str10-deep
   2    BL10M    buried-locking-str10-middle
   3    BL10S    buried-locking-str10-shallow
   4    BL13D    buried-locking-mur13-deep

"""

from numpy import *
import os,sys,glob

from clawpack.clawutil.util import fullpath_import
from clawpack.clawutil import multip_tools, clawmultip_tools

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

common_code_dir = os.path.join(CHT, 'common_code')
cases_dtopos = fullpath_import(f'{common_code_dir}/cases_dtopos.py')

dry_run = True  # If True, only print out settings, do not run GeoClaw
#dry_run = False  # If True, only print out settings, do not run GeoClaw

# what to do:
run_code = True
make_plots = True


# location for big files for different computer environments:
this_dir = os.getcwd()
HOME = os.environ['HOME']

if 'rjl/git' in this_dir:
    computer = 'rjl-laptop'
    scratch_dir = this_dir.replace('rjl/git', 'rjl/scratch')

elif '/mmfs1/home' in this_dir:
    computer = 'hyak'
    scratch_dir = this_dir.replace('/mmfs1/home', '/gscratch/tsunami')

elif '/home1' in this_dir:
    computer = 'tacc'
    #scratch_dir = this_dir.replace('/home1', '/scratch')
    try:
        SCRATCH = os.environ['SCRATCH']
        scratch_dir = this_dir.replace(HOME, SCRATCH)
    except:
        scratch_dir = this_dir  # if $SCRATCH not set

else:
    computer = 'unknown'
    scratch_dir = this_dir

# where to find all the dtopo files: (modified below on TACC)
if 1:
    # original kinematic ruptures from SPECFEM3D:
    instant = False
    dtopo_dir = f'{CHT}/dtopo/CSZ_groundmotions/dtopo30sec/dtopofiles'
else:
    # static ruptures without subevents, using Okada:
    instant = True
    dtopo_dir = f'{CHT}/dtopo/CSZ_groundmotions/nosubevents_251229/dtopofiles'


if computer == 'tacc':
    #dtopo_dir = dtopo_dir.replace('/home1', '/scratch')
    dtopo_dir = dtopo_dir.replace(HOME, '/work2/04137/rjl/CHTshare/')


print('scratch_dir = ',scratch_dir)
print('dtopo_dir = ',dtopo_dir)

# where to put output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)

# create scratch directory for output, if it doesn't exist:
os.system('mkdir -p %s' % runs_dir)

# path to geoclaw executable:
# (should agree with how EXE is set in Makefile used to compile)
if run_code:
    xgeoclaw_path = f'{CHT}/geoclaw_runs/xgeoclaw-v5.13.1'
    if computer == 'tacc':
        xgeoclaw_path = '/work2/04137/rjl/CHTshare/clawpack-share/tacc/xgeoclaw_251229'
else:
    xgeoclaw_path = None  # do not run GeoClaw code


# Specify the list of dtopo files to loop over for geoclaw runs:

# List of all events from CoPes Hub ground motions:
# For naming and numbering convention, see
#   https://depts.washington.edu/ptha/CHTuser/dtopo/groundmotions/
# Note that there may be different versions of these events, so which
# version is used depends on what directory dtopo_dir points to

try:
    import CHTtools
    all_events = CHTtools.all_events()  # returns the list of 36 events
except:
    # if CHTtools isn't found, construct the list of 36 events here:

    depths = ['D','M','S']

    # buried_locking events:
    all_events = [f'BL10{depth}' for depth in depths] \
               + [f'BL13{depth}' for depth in depths] \
               + [f'BL16{depth}' for depth in depths] \

    all_events += [e.replace('L','R') for e in all_events]  # add random events
    all_events += [e.replace('B','F') for e in all_events]  # add ft events

    all_events.sort()

# select which subset of all_events to run at command line
# or in slurm script by specifying NPROCS FIRST_EVENT LAST_EVENT

# if dtopo_dir points to a directory that has instantaneous versions
# (static displacement rather than kinematic time-dependent)
# then you could set `instant = True` to use these, provided they
# have file names such as BL10D_instant.dtt3 (with the same numbering 1-36):

if instant:
    all_events = [e+'_instant' for e in all_events]



if __name__ == '__main__':

    import sys

    # see the doc string at the top of this file for more details on
    # how to run this script

    # parse command line arguments:

    # always expect at least one argument NPROCS, how many to run in parallel
    try:
        nprocs = int(sys.argv[1])
    except:
        raise Exception('*** Missing integer argument nprocs on command line')

    # might also have arguments FIRST_EVENT LAST_EVENT
    # in which case restrict all_events set above to these events
    # (with e.g. 1 4 meaning the first 4 events, i.e. all_events[0:4]):

    if len(sys.argv) == 3:
        msg = '*** unexpected number of command line arguments'
        raise ValueError(msg)

    if len(sys.argv) == 4:
        # user wants a subset of all_events:
        ievent_first = int(sys.argv[2]) - 1
        ievent_last = int(sys.argv[3])
        if 0 <= ievent_first <= ievent_last <= len(all_events):
            events = all_events[ievent_first:ievent_last]
        else:
            msg = f'*** expected 0 < FIRST_EVENT={sys.argv[2]}' \
                    +f' <= LAST_EVENT={sys.argv[3]}' \
                    +f' <= {len(all_events)} = len(all_events) in arguments'
            raise ValueError(msg)
    else:
        # FIRST_EVENT LAST_EVENT not specified on command line:
        events = all_events

    dtopo_files = ['%s/%s.dtt3' % (dtopo_dir,e) for e in events]

    print('\n--------------------------')
    if dry_run:
        # just print out settings, no runs...
        print('DRY RUN - settings in runclaw_makeplots_dtopos.py')

    print('Will run GeoClaw for %i dtopo files' % len(dtopo_files))
    print('dtopo files should be in dtopo_dir:\n    ', dtopo_dir)
    print('list of dtopo_files to process:')
    #for fname,fpath in zip(dtopo_names, dtopo_files):
    for fpath in dtopo_files:
        print('  %s' %  os.path.split(fpath)[-1])
        if not os.path.isfile(fpath):
            print('    *** file not found')
    print('output will go in \n    %s/geoclaw_outputs/' % runs_dir)
    print('plots will go in \n    %s/geoclaw_plots/' % runs_dir)
    print('nprocs = %i jobs will run simultaneously' % nprocs)
    print('OMP_NUM_THREADS = ', os.environ['OMP_NUM_THREADS'])
    print('GeoClaw executable:\n    ',xgeoclaw_path)

    if dry_run:
        print('Set dry_run=False and re-execute to run GeoClaw')
        print('--------------------------')
    else:
        foutfile = f'{runs_dir}/geoclaw_outputs/fortran_output.txt'
        poutfile = f'{runs_dir}/geoclaw_outputs/python_output.txt'
        print(f'Fortran output will go to \n    {foutfile}')
        print(f'Python output will go to \n    {poutfile}')

        # make list of dictionaries with parameters for each case:
        caselist = cases_dtopos.make_all_cases_dtopos(dtopo_dir, dtopo_files,
                                        runs_dir, xgeoclaw_path, make_plots)

        # run all cases using nprocs processors:
        multip_tools.run_many_cases_pool(caselist, nprocs,
                                         clawmultip_tools.run_one_case_clawpack)
