#!/bin/bash
#SBATCH -J job-run          # Job name
#SBATCH -o job-run.o%j      # Name of stdout output file
##SBATCH -e job-run.e%j      # Name of stderr error file
#SBATCH -p skx-dev          # Queue (partition) name
##SBATCH -p skx              # Queue (partition) name
#SBATCH -N 1                # Total # of nodes 
#SBATCH -n 1                # Total # of mpi tasks
#SBATCH -t 2:00:00         # Run time (hh:mm:ss)
##SBATCH -t 15:00:00         # Run time (hh:mm:ss)
##SBATCH --mail-user=username@tacc.utexas.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A DS-portal-lma3   # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

date # record start time


export OMP_NUM_THREADS=48  # for a single run with 1 dtopo file

# Note: executable xgeoclaw to use for Fortran is set by Makefile
# These are not needed, since not compiling fortran:
#export FC=ifx
#export FFLAGS='-O2 -qopenmp'

export CHT=$HOME/CopesHubTsunamis  # for this project

# version of geoclaw to use for Python code:
export CLAW=/work2/04137/rjl/CHTshare/clawpack-share

export PYTHONPATH=$CLAW

date # record beginning time

pwd
module list
echo CHT=$CHT
echo CLAW=$CLAW
echo PYTHONPATH=$PYTHONPATH

# commands to run:
#make scrdir  # create scratch directory
#make data    # using setrun_case.py
#make output  # using EXE set in Makefile
make plots   # using setplot_case.py

date # record end time
