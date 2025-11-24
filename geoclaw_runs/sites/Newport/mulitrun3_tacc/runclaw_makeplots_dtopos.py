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

from clawpack.clawutil.util import fullpath_import
from clawpack.clawutil import multip_tools, clawmultip_tools

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

common_code_dir = os.path.join(CHT, 'common_code')
cases_dtopos = fullpath_import(f'{common_code_dir}/cases_dtopos.py')
plot_gauges_site = fullpath_import(f'{common_code_dir}/plot_gauges_site.py')

if 0:
    # can now import from clawutil?
    clawmultip_tools = fullpath_import(f'{common_code_dir}/clawmultip_tools.py')
    multip_tools = fullpath_import(f'{common_code_dir}/multip_tools.py')

#dry_run = True  # If True, only print out settings, do not run GeoClaw
dry_run = False  # If True, only print out settings, do not run GeoClaw

# what to do:
run_code = True
make_plots = False


# location for big files for different computer environments:
this_dir = os.getcwd()

if 'rjl/git' in this_dir:
    computer = 'rjl-laptop'
    scratch_dir = this_dir.replace('rjl/git/CopesHubTsunamis/geoclaw_runs', \
                                   'rjl/scratch/CHT_runs')

elif '/mmfs1/home' in this_dir:
    computer = 'hyak'
    scratch_dir = this_dir.replace('/mmfs1/home', '/gscratch/tsunami')

elif '/home1' in this_dir:
    computer = 'tacc'
    scratch_dir = this_dir.replace('/home1', '/scratch')

else:
    computer = 'unknown'
    scratch_dir = this_dir

# where to find all the dtopo files:
dtopo_dir = f'{CHT}/dtopo/CSZ_groundmotions/dtopo30sec/dtopofiles'

if computer == 'tacc':
    dtopo_dir = dtopo_dir.replace('/home1', '/scratch')


print('dtopo_dir = ',dtopo_dir)

# where to put output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)

# create scratch directory for output, if it doesn't exist:
os.system('mkdir -p %s' % runs_dir)

# path to geoclaw executable:
# (should agree with how EXE is set in Makefile used to compile)
if run_code:
    xgeoclaw_path = f'{CHT}/geoclaw_runs/xgeoclaw_v5-13-1'
    if computer == 'tacc':
        xgeoclaw_path = f'{CHT}/geoclaw_runs/xgeoclaw_v5-13-1_ifx'
else:
    xgeoclaw_path = None  # do not run GeoClaw code

# number of events to run and/or plot simultaneously:
#nprocs = 9 # now set on command line (or in slurm script)

# Specify the list of dtopo files to loop over for geoclaw runs:

if 0:
    # all files found
    all_dtopo_files = glob.glob('%s/*tt3' % dtopo_dir)
    dtopo_files = all_dtopo_files
    # not the instant deformation versions:
    dtopo_files = [f for f in all_dtopo_files if 'instant' not in f]
    dtopo_files.sort()

# specify events...

depths = ['D','M','S']

# buried_locking events:
all_events = [f'BL10{depth}' for depth in depths] \
           + [f'BL13{depth}' for depth in depths] \
           + [f'BL16{depth}' for depth in depths] \

all_events += [e.replace('L','R') for e in all_events]  # add random events
all_events += [e.replace('B','F') for e in all_events]  # add ft events

events = all_events
events.sort()

#events = events[:9]
#events = events[9:]
#events = events[4:12]

events = ['BL10D'] + events[12:15]

instant = False
if instant:
    events = [e+'_instant' for e in events]


dtopo_files = ['%s/%s.dtt3' % (dtopo_dir,f) for f in events]


dtopo_names = []
for f in dtopo_files:
    dtopo_name = os.path.splitext(os.path.split(f)[-1])[0]
    dtopo_names.append(dtopo_name)

#print('dtopo_files = ',dtopo_files)
#print('dtopo_names = ',dtopo_names)

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
    #for fname,fpath in zip(dtopo_names, dtopo_files):
    for fpath in dtopo_files:
        print('  %s' %  os.path.split(fpath)[-1])
        if not os.path.isfile(fpath):
            print('    *** file not found')
    print('output will go in \n    %s/geoclaw_outputs/' % runs_dir)
    print('nprocs = %i jobs will run simultaneously' % nprocs)
    print('OMP_NUM_THREADS = ', os.environ['OMP_NUM_THREADS'])
    print('GeoClaw executable:\n    ',xgeoclaw_path)

    if dry_run:
        print('Set dry_run=False and re-execute to run GeoClaw')
        print('--------------------------')
    else:
        # make list of dictionaries with parameters for each case:
        caselist = cases_dtopos.make_all_cases_dtopos(dtopo_dir, dtopo_files,
                                        runs_dir, xgeoclaw_path, make_plots)

        # run all cases using nprocs processors:
        multip_tools.run_many_cases_pool(caselist, nprocs,
                                         clawmultip_tools.run_one_case_clawpack)
