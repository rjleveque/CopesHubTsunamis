#!/bin/bash
#SBATCH -J job-test           # Job name
#SBATCH -o job-test.o%j       # Name of stdout output file
#SBATCH -e job-test.e%j       # Name of stderr error file
#SBATCH -p skx-dev         # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-rjl    # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

date # record start time

# not needed, since not compiling fortran:
#export FC=ifx
#export FFLAGS='-O2 -qopenmp'

export CHT=$HOME/CopesHubTsunamis  # for this project

# set CLAW used for python
export CLAW=/work2/04137/rjl/CHTshare/clawpack-share
export PYTHONPATH=$CLAW

pwd
module list
echo CHT=$CHT
echo CLAW=$CLAW
echo PYTHONPATH=$PYTHONPATH

# commands to run:
python make_allevents_kmz.py

date # record end time
