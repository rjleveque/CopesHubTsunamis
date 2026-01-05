
Create the executable xgeo_computeB in this directory $CHT/common_code via:

    make -f Makefile_computeB

Then in any directory with a setrun.py that defines gauges:

    make data  # creates gauges.data, topo.data, etc 
    $CHT/common_code/xgeo_computeB
   
This will print out something like this for each gauge:

 ==============================================================================
Gauge    103  At x=-124.1704167        y=  47.1162500
          xcenter =-124.1704167 ycenter =  47.1162500 B in cell =    -0.431
 
                |                    |                   |               |
  47.1163889  --+--------------------+-------------------+---------------+--
                |       -0.443       |       -0.432      |       -0.421  |
  47.1162963  --+--------------------+-------------------+---------------+--
                |       -0.413       |       -0.431      |       -0.428  |
  47.1162037  --+--------------------+-------------------+---------------+--
                |       -0.416       |       -0.439      |       -0.442  |
  47.1161111  --+--------------------+-------------------+---------------+--
                |                    |                   |               |
         -124.1705556           -124.1704630      -124.1703704     -124.1702778
 ==============================================================================

The first line gives the gauge location as read from gauges.data
The next line gives the cell center for the grid cell this gauge lies in
(currently hard-wired to assume dx=dy=1/3" resolution, starting from xlower,ylower
as read from claw.data (which was specified in setrun.py)

The grid then shows the GeoClaw cell-averaged topo values B(i,j) in a 3x3
array of grid cells (at resolution dx,dy), as computed by the same method as
used in a GeoClaw run (i.e. all the topofiles specified in topo.data are read in
and a piecewise bilinear function computed from these, followed by exact
integration of this function over each grid cell).

The center grid cell in this 3x3 grid contains the gauge.

Currently it prints the cell-averaged B value in each cell.
In GeoClaw runs, the gauge values printed at each time are either:
    - Interpolated between 4 cell centers surrounding the gauge if all 4 have h>0
    - The cell-averaged values of h,hu,hv,eta=h+B if any of these 4 cells is dry
If the gauge is at the cell center, both approaches give the same value.

Still to do:

Accept a sea_level parameter and compute the interpolate B value properly
in the case where all 4 surrounding cell centers have B < sea_level (i.e.
initially all are wet with this value of sea_level).
(And/or use sea_level as set in setrun.py)

Accept dx,dy values so the topography can be examined at different resolutions.
The finest resolution for the run could be computed from info in setrun.py
(the coarse level resolution, amr_levels_max, and refinement_ratios), but some
gauges may not be on finest-level grids, or we may want to easily see how changing
the resoluton affects things.

     

