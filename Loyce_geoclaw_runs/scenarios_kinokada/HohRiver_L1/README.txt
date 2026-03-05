Set up to run a single dtopo file on TACC

Submit batch job via:
    sbatch run_1dtopo.slurm
This script defines the version of CLAW being used
Uses 48 OpenMP threads for a single run
Check which queue is used and time limit for job run

This uses Makefile, which sets
    EXE to executable in $CHTshare
    SCRATCH to the same path as this directory, with /home1 replaced by /scratch
    OUTDIR to $SCRATCH/_output
    PLOTDIR to $SCRATCH/_plots

Note: it is NOT currently set up to make directories
    .../geoclaw_outputs/dtopo_name
    .../geoclaw_plots/dtopo_name

setrun_case.py sets the case parameters in the "main" program at the bottom

setplot_case.py sets the case parameters near the top before defining the
     setplot() function

