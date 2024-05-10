"""
Code to perform multiple GeoClaw runs, with a different dtopo file
for each.

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

sys.path.insert(0, os.path.join(root_dir, 'common_code'))
from multip_tools import make_all_cases_dtopos, run_one_case_dtopo
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
dtopo_dir = os.path.join(root_dir, 'dtopo/CSZ_groundmotions/dtopofiles')
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
xgeoclaw_path = os.path.join(root_dir, 'geoclaw_runs/xgeoclaw')

# number of GeoClaw jobs to run simultaneously:
nprocs = 6

# Specify the list of dtopo files to loop over for geoclaw runs:

if 1:
    # all files found
    all_dtopo_files = glob.glob('%s/*tt3' % dtopo_dir)
    # not the instant deformation versions:
    dtopo_files = [f for f in all_dtopo_files if 'instant' not in f]
    dtopo_files.sort()

if 1:
    # test on only a few files from list above:
    dtopo_files = dtopo_files[12:18] + dtopo_files[0:6]
    
if 0:
    # or set a list of dtopo_names by other means, e.g.
    events = ['buried-locking-mur13-deep','buried-locking-str10-shallow']
    dtopo_files = ['%s/%s.dtt3' % (dtopo_dir,f) for f in events]

dtopo_names = []
for f in dtopo_files:
    dtopo_name = os.path.splitext(os.path.split(f)[-1])[0]
    dtopo_names.append(dtopo_name)

#print('dtopo_files = ',dtopo_files)
#print('dtopo_names = ',dtopo_names)

if __name__ == '__main__':
    
    print('\n--------------------------')
    if dry_run:
        # just print out settings, no runs...
        print('DRY RUN - settings in run_geoclaw_dtopos.py')

    print('Will run GeoClaw for %i dtopo files' % len(dtopo_files))
    print('list of dtopo_files to process:\n', dtopo_files)
    print('dtopo files should be in dtopo_dir:\n    ', dtopo_dir)
    print('output will go in \n    %s/geoclaw_outputs/' % runs_dir)
    print('nprocs = %i jobs will run simultaneously' % nprocs)
    print('GeoClaw executable:\n    ',xgeoclaw_path)

    if dry_run:
        print('Set dry_run=False and re-execute to run GeoClaw')
        print('--------------------------')
    else:    
        # make list of dictionaries with parameters for each case:
        caselist = make_all_cases_dtopos(dtopo_dir, dtopo_files, xgeoclaw_path,
        				 runs_dir)

        # run all cases using nprocs processors:
        run_many_cases_pool(caselist, nprocs, run_one_case_dtopo)
