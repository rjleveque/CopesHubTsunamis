#!/bin/bash
#SBATCH -J job-test           # Job name
#SBATCH -o job-test.o%j       # Name of stdout output file
#SBATCH -e job-test.e%j       # Name of stderr error file
#SBATCH -p skx-dev         # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 00:10:00        # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-rjl    # Allocation name (req'd if you have more than 1)

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

pwd
module list
echo CHT=$CHT
echo CLAW=$CLAW
echo PYTHONPATH=$PYTHONPATH

# commands to run:
python runclaw_makeplots_dtopos.py $NPROCS

date # record end time
