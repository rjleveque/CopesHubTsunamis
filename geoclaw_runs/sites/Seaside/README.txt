
# Seaside, OR modeling

Set environment variable CHT to top level of this CopesHubTsunami repository.

Running the code requires topo and dtopo files.  TODO: add instructions

Note that some topo or dtopo files use longitude E and need to be shifted by
subtracting 360 from the xlower value in the header.  For this problem the
topo files are fine but the dtopo file

    $COPES/dtopo/dtopofiles/CSZ_SM1.tt3
    
must be modified by changing the fourth line from

    2.33000000000000e+02   xlower
    
to

    -1.27000000000000e+02   xlower
    
Otherwise the seafloor deformation will be outside of the computational
domain and no tsunami will be generated.

Then

    cd $CHT/geoclaw_runs/sites/seaside/CSZ_SM1
    make data

will produce .data files and also .kml files to view on Google Earth to see
the extents of topo and dtopo files used, AMR flagregions, gauge locations, etc.

To run the code:

    make .output

--------------------------------------------------------
Making plots:

The Python scripts in `$CHT/geoclaw_runs/sites/seaside` that are described
below will create plots and an animation in

    $CHT/geoclaw_runs/sites/seaside/CSZ_SM1/_plots
    
or more generally in a _plots subdirectory within the run directory, if
used on a different event than CSZ_SM1.
    
To make gauge plots:

    cd $CHT/geoclaw_runs/sites/seaside/CSZ_SM1
    python ../plot_gauges.py

To make an animation using data in the fgout files:

    python ../make_fgout_animation.py

--------------------------------------------------------
To make fgmax plots and kmz file:

It is necessary to first insure that you have a file
    $CHT/geoclaw_runs/sites/seaside/fgmax0001_13s_B0.asc
that has the pre-seismic topography B0 at all fgmax points.
This file is created by running the code in the seaside1_B0 subdirectory:

    cd $CHT/geoclaw_runs/sites/seaside/seaside1_B0/
    make .output
    
This directory contains a setrun.py file that uses the same refinement levels
and regions as in CSZ_SM1, but with no dtopo file specified so there is no
uplift or subsidence.  It is also set to run for only a few seconds and
produce no time frame or gauge output, only _output/fgmax0001.txt
This captures the GeoClaw topography values B on the fgmax grid.  This file
is converted to the form of a topofile via:

    cd $CHT/geoclaw_runs/sites/seaside
    python make_B0.py
    
which produces
    $CHT/geoclaw_runs/sites/seaside/fgmax0001_13s_B0.asc

Note that if the fgmax grid is changed, or the topofiles, refinement levels,
regions, etc. car changed, a new B0 file will need to be created.
Also note that it is important to check that all the refinement regions
around the fgmax grid are turned on to their maximum level starting at t=0
in order to capture B0 at the same resolution as the fgmax values computed in
a job run for some event.
    
Once the B0 file exists, you can create fgmax plots via:

    cd $CHT/geoclaw_runs/sites/seaside/CSZ_SM1
    python ../process_fgmax.py
    

--------------------------------------------------------
Adapting to a different site:
    
The scripts can be modified by changing `event` to be
used for additional runs using different dtopo files.

To adapt to a different location, copy the script, change `location`,
and you may also need to make some other modifications to better
suit the new location. (e.g. figsize based on the aspect ratio of the fgmax
plot.

You will also have to redefine what points to consider `onshore`,
since for Seaside we set it to consider any point with x > -123.93 onshore
in order to capture the harbor and rivers also in the plot of "depth".

What we really plot for "depth" is what we call zeta, equal to:

    zeta = h "onshore"
         = h + B0 "offshore"
         
Normally we define

    onshore = where B0 > 0
    
This is all based on the pre-seismic topography B0 (elevation relative to MHW).

We always define the water elevation (relative to MHW) at any time before or
during the tsunami as

    eta = h + B
    
as is standard in tsunami models, and is based on the elevation B at that time
(which possibly includes subsidence or uplift after the earthquake happens),
but using h+B0 for the plots together with h where the land was originally dry
has the advantage that zeta is continuous at the
original shoreline (where B0=0).

zeta has the interpretation of the apparent water elevation in a harbor
or river to an observer standing on shore, since the observer would drop
with co-seismic subsidence along with the land and the initial water
surface before the tsunami arrives. So, for example, if the tsunami
raises the water level in a harbor by 1 meter relative to the docks at
some time (i.e. the depth h has increased by 1 meter at a point where
initially h0+B0=0), then zeta would be 1 meter, as makes sense to an
onshore observer, whereas if the harbor and land all dropped by
B - B0 = -3 meters of subsidence, then eta = h+B = -2 m in the harbor,
which is less useful to know.
