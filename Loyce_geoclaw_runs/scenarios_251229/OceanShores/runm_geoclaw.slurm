#!/bin/bash
#SBATCH -J job-run          # Job name
#SBATCH -o job-run.o%j      # Name of stdout output file
#SBATCH -e job-run.e%j      # Name of stderr error file
##SBATCH -p skx-dev          # Queue (partition) name
#SBATCH -p skx               # Queue (partition) name
#SBATCH -N 1                # Total # of nodes 
#SBATCH -n 1                # Total # of mpi tasks
#SBATCH -t 6:00:00         # Run time (hh:mm:ss)
##SBATCH -t 00:10:00         # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-lma3   # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

date # record start time

# Set NPROCS from sbatch command and compute OMP_NUM_THREADS
# so that NPROCS * OMP_NUM_THREADS <= 48 (number of cores per node)

export NPROCS=$1
if [ "$NPROCS" = "" ]; then
    echo "SLURM ERROR -- Must specify NPROCS in sbatch command"
    exit 1
elif [ "$NPROCS" -gt 48 ]; then
    echo "SLURM ERROR -- Must specify NPROCS <= 48"
    exit 1
else
    export OMP_NUM_THREADS=$(( 48 / NPROCS ))
    echo NPROCS = $NPROCS
    echo OMP_NUM_THREADS = $OMP_NUM_THREADS
fi

export FIRST_EVENT=$2
export LAST_EVENT=$3

# not needed, since not compiling fortran:
#export FC=ifx
#export FFLAGS='-O2 -qopenmp'

export CHT=$HOME/CopesHubTsunamis  # for this project

# Note: executable xgeoclaw to use for Fortran is set by
#       xgeoclaw_path in runclaw_makeplots_dtopos.py

# version of geoclaw to use for Python code:
#export CLAW=$HOME/clawpack_src/clawpack-v5.13.1
export CLAW=/work2/04137/rjl/CHTshare/clawpack-share

export PYTHONPATH=$CLAW

date # record beginning time

pwd
module list
echo CHT=$CHT
echo CLAW=$CLAW
echo PYTHONPATH=$PYTHONPATH

# commands to run:
python runclaw_makeplots_dtopos.py $NPROCS $FIRST_EVENT $LAST_EVENT

date # record end time
