
Sample code to test running GeoClaw on TACC, looping over multiple dtopos

Files to modify for a set of runs:


----------------------------------------
setrun_case.py
----------------------------------------

This file contains the setrun function that will be run for each dtopo
to create the *.data files needed by the Fortran code.

Note that this setrun function has an additional argument `case`
that is not usually in setrun.  This is a dictionary that will be passed
in when creating the `*.data` files for a particular case.
For our situation, the only important entry in the dictionary is
`dtopofiles` since the value `case['dtopofiles']` will be used in setting
`dtopo_data.dtopofiles`.  (Also `case['outdir']` is used if doing a
restart.)

Everything else in `setrun` will be exactly the same for all events being run.

Note that when looping over many dtopo events, you may not want to create
time frame output of the full AMR solution at lots of output times, and
perhaps not at any output times, since this output tends to be huge and may
not be very useful for post-processing the results.  Instead one often only
wants one or more of the following:
 - gauge output: time series at specified synthetic gauge locations,
 - fgmax output: maximum inundation depth, speed, etc. on a fixed grid,
   usually covering a small portion of the computational domain,
 - fgout output: snapshots of the solution on fixed uniform grid, possibly
   at many times, but only over a small portion of the computational domain.

The example included here is set up to produce no time frame output and only
record time series at a single gauge (over a very short run for illustration).

----------------------------------------
setplot_case.py
----------------------------------------

This file is only needed if `make_plots == True` in 
`runclaw_makeplots_dtopos.py` (see below).

This file can then be similar to a standard `setplot.py` that specifies how
the plots should be made for each time frame of output data. However,
as noted above, you may not want to save full AMR time frame results.
So alternatively (or in addition), the setplot function can call any other
post-processing scripts that you want to call after performing the geoclaw
runs, for example to make special plots for each gauge, to plot fgmax
results, or to make animations from fgout grid data.

Similar to `setrun_case`, the setplot function in this file has an
additional argument `case` that is used in the post-processing.

The example included here creates one gauge plot and a text file in the
`_plots` directory for each event.

----------------------------------------
runclaw_makeplots_dtopos.py
----------------------------------------

set dry_run = True to just print out info about what will be done   
              False to actually run the code and/or make plots

Check that dtopo_dir and topo_dir will be set properly

Output and plots will be sent to scratch_dir
which is constructed to be the full path within CopesHubTsunamis but in
$SCRATCH instead of in $HOME.

For example, if the directory containing this code is

    /home1/04137/rjl/CopesHubTsunamis/geoclaw_runs/tacc-test

then the output/plots will be sent to

    /scratch/04137/rjl/CopesHubTsunamis/geoclaw_runs/tacc-test/geoclaw_outputs
    /scratch/04137/rjl/CopesHubTsunamis/geoclaw_runs/tacc-test/geoclaw_plots

within these directories, there will be one subdirectory for each event,
with names like `_output_BL10D` or `_plots_BL10D`.

Set `run_code = True` to run geoclaw for each event
Set `make_plots = True` to make plots (or other post-processing) for each event

If you have already run geoclaw on each events and only want to redo the
post-processing, make sure you set `run_code = False`.

If `run_code` is True, make sure `xgeoclaw_path` is properly set to the
executable (which must have been previously compiled using an appropriate
Makefile and version of GeoClaw).

Set events to the list of event names to be run.
Event names like 'BL10D' are now typically expected, and the dtopofile should be
in f'{dtopo_dir}/{event}.dtt3'

The sample code sets all_events to a list of 36 events and then selects only
the first 2 for this test case.

After modifying this file, test it by setting `dry_run = True` and then:

    python runclaw_makeplots_dtopos.py 2

The number 2 indicates that it should do 2 runs at a time using Python
multiprocessing tools.

If the screen output from this looks ok, change to `dry_run = False`
and submit a batch run using slurm.

----------------------------------------
slurm script for batch submission
----------------------------------------

The script `$CHT/tacc_stampede3/runm_geoclaw-test.slurm` should be
modified to specify your allocation in place of:

    #SBATCH -A DS-portal-rjl   # Allocation name (req'd if you have more than 1)

This script is set up to run 2 jobs in parallel.  Since each node on
stampede3 has 48 cores, it is set to use `OMP_NUM_THREADS=24` for each
job.  Other scripts in `$CHT/tacc_stampede3/` are set up to run e.g.
6 jobs in parallel with 8 threads for each.

Other things you might want to change:

    #SBATCH -p skx-dev         # Queue (partition) name

The `skx-dev` queue is for short test jobs, for long runs you will need to
use `skx` instead, and increase the allowed run time:

    #SBATCH -t 00:10:00        # Run time (hh:mm:ss)


Submit a batch job using the `sbatch` command.




