#!/bin/bash
#SBATCH -J job-okada           # Job name
#SBATCH -o job-okada.o%j       # Name of stdout output file
#SBATCH -e job-okada.e%j       # Name of stderr error file
#SBATCH -p skx             # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 05:00:00        # Run time (hh:mm:ss)  # max 2hrs for skx-dev
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-rjl    # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

#export FC=gfortran
#export FFLAGS='-O2 -fopenmp'
#export OMP_NUM_THREADS=6
#export NPROCS=8  # NPROCS * OMP_NUM_THREADS = 48
#export CLAW=/home1/04137/rjl/clawpack_src/clawpack-v5.12.0
#export PYTHONPATH=$PYTHONPATH:$CLAW
export CHT=$HOME/CopesHubTsunamis  # for this project
#export CHTuser=/corral/projects/NHERI/projects/7f2e74be-d7ca-4e0e-b69a-22c24840b078/CHTuser

date
pwd
module list
echo CHT=$CHT

python OkadaFromSubsampledSlip.py

date
