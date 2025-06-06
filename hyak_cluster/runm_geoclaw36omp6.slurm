#!/bin/bash

# JOB NAME
#SBATCH --job-name=runm36omp6

# ALLOCATION DEFINITION
# The account and partition options should be the same
# except in a few cases (e.g. ckpt queue, genpool queue)
#SBATCH --account=tsunami
#SBATCH --partition=compute

# RESOURCES
#SBATCH --nodes=1  # total number of nodes allocated
#SBATCH --ntasks-per-node=36  # cores per node

# WALL TIME
# Do not specify a wall time significantly more than your job needs
# Common acceptable time formats:
#    hours:minutes:seconds e.g. 3:00:00 for 3 hours
#    minutes
#    days-hours
#SBATCH --time=10:00:00

# MEMORY PER NODE
# See above "Specifying Memory Size" for options
#SBATCH --mem=100G  # e.g. --mem=100G for 100 GB of memory

# WORKING DIRECTORY ENTRYPOINT
# Specify the working directory for this job
#SBATCH --output=job-runm36.%j.out
#SBATCH --error=job-runm36.%j.err

# Turn on email notifications
##SBATCH --mail-type=ALL
##SBATCH --mail-user=rjl@uw.edu

# Run the commands to run your program here
# e.g. load modules, copy input.output files, run program, etc.

#module avail
module load intel/oneAPI/2021.1.1

# Export all environment variables to the batch job session
#SBATCH --export=all

# Set or reset any environment variables here:
source $HOME/.bashrc   # needed for conda setup

export FC=ifort
export FFLAGS='-O2 -qopenmp'
export OMP_NUM_THREADS=6
export CLAW=/gscratch/tsunami/clawpack-v5.12.0
export PYTHONPATH=$PYTHONPATH:$CLAW
export CHT=$HOME/CopesHubTsunamis  # for this project

date
pwd
echo CHT=$CHT
conda activate /gscratch/tsunami/miniconda3/envs/geo5
#python run_geoclaw_dtopos.py
python runclaw_makeplots_dtopos.py
date

