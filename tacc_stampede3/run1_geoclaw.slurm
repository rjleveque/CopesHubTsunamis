#!/bin/bash
#SBATCH -J job-run1           # Job name
#SBATCH -o job-run1.o%j       # Name of stdout output file
#SBATCH -e job-run1.e%j       # Name of stderr error file
#SBATCH -p skx             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 05:00:00        # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-rjl    # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

date # record start time


export FC=ifx
export FFLAGS='-O2 -qopenmp'

export CHT=$HOME/CopesHubTsunamis  # for this project

export CLAW=/work2/04137/rjl/CHTshare/clawpack-share
export PYTHONPATH=$CLAW
export OMP_NUM_THREADS=48

pwd
module list
echo CHT=$CHT
echo CLAW=$CLAW
echo PYTHONPATH=$PYTHONPATH

# commands to run:
make scrdir
make data
make output

date # record end time
