"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import
from __future__ import print_function
import os
import numpy as np
from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgout_tools, fgmax_tools



try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

root_dir = os.environ['CHT']
input_files_dir = root_dir + '/input_files/'
print('setting input_files_dir = ',input_files_dir)

rundir = os.getcwd()
print('rundir = %s' % rundir)

if '/projects' in rundir:
    topodir = '/projects/rale6846/topo/topofiles'  # on CU
    dtopodir = '/projects/rale6846/dtopo/dtopofiles'  # on CU
else:
    topodir = '/Users/rjl/topo/topofiles'        # on Randys laptop
    #dtopodir = '/Users/rjl/B/dtopo/dtopofiles'   # on Randys laptop
    dtopodir = root_dir + '/dtopo/CSZ_groundmotions'   # on Randys laptop


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:


    clawdata.lower[0] = -135.      # west longitude
    clawdata.upper[0] = -122.      # east longitude
    clawdata.lower[1] = 38.5       # south latitude
    clawdata.upper[1] = 53.5       # north latitude

    clawdata.num_cells[0] = 13
    clawdata.num_cells[1] = 15


    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False   # True to restart from prior results
    clawdata.restart_file = ''

    tstart_finestgrid = 0. #14*60.

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        ## ADJUST:
        clawdata.num_output_times = 4
        clawdata.tfinal = 2*3600.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = []

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True
        

    clawdata.output_format = 'binary'

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # negative checkpoint_style means alternate between aaaaa and bbbbb files
    # so that at most 2 checkpoint files exist at any time, useful when
    # doing frequent checkpoints of large problems.

    clawdata.checkpt_style = -2

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = 3600.*np.arange(1,16,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 1deg, 6', 3', 45", 9", 3", 1", 1/3"":
    # refinement_ratios = [10,2,4,5,3,3,3]

    # dx = dy = 1deg, 6', 1', 15"
    refinement_ratios = [10,6,4]
    amrdata.refinement_ratios_x = refinement_ratios
    amrdata.refinement_ratios_y = refinement_ratios
    amrdata.refinement_ratios_t = refinement_ratios



    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.5

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    #topodir = '/Users/rjl/topo/topofiles/'
    topofiles = topo_data.topofiles

    # 1-minute topo:
    topofiles.append([3, topodir + '/etopo1_-163_-122_38_63.asc'])

    if 0:
        # 2-second topo:
        topofiles.append([3, topodir + '/PT_2sec_center.asc'])
        topofiles.append([3, topodir + '/PS_2sec_center.asc'])
        topofiles.append([3, topodir + '/SJdF_2sec_center.asc'])

    if 0:
        # 1/3-second topo:
        topofiles.append([3, topodir + '/GH_13sec.asc'])
        topofiles.append([3, topodir + '/WB_13sec.asc'])

    #topofiles.append([3, '/Users/rjl/topo/WA/astoria_13_mhw_2012/GH_13sec.asc'])
    #topofiles.append([3, '/Users/rjl/topo/WA/astoria_13_mhw_2012/WB_13sec.asc'])


    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]

    dtopo_data.dtopofiles = [[3, dtopodir + '/buried-locking-skl16-shallow.dtt3']]

    dtopo_data.dt_max_dtopo = 10


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    rundata.qinit_data.variable_eta_init = True  # for subsidence


    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    flagregions = rundata.flagregiondata.flagregions  # initialized to []
    
    #RRdir = '/Users/rjl/git/paleo2020/topo/RuledRectangles'
    RRdir = root_dir + '/topo/regions/'
    
    # Computational domain Variable Region - 1degree to 6min to 3min:
    # Level 3 below is 3 min
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [clawdata.lower[0]-0.1,
                                 clawdata.upper[0]+0.1,
                                 clawdata.lower[1]-0.1,
                                 clawdata.upper[1]+0.1]
    flagregions.append(flagregion)

    # Source Variable Region - 3min to 45sec:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_dtopo'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 30.
    flagregion.spatial_region_type = 1  # Rectangle
    source_region = [-129,-122.1,38.5,50]
    flagregion.spatial_region = source_region
    flagregions.append(flagregion)

    # Continential shelf Variable Region - 45sec:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.*3600.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51.data')
    flagregions.append(flagregion)
    
    # Continential shelf Variable Region - 45sec:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.*3600.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46.data')
    flagregions.append(flagregion)

    if 0:
    
        # Source Variable Region - 45sec:
        # to fill in coast missing from ruled rectangles above
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_north_coast'
        flagregion.minlevel = 4
        flagregion.maxlevel = 4
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        source_region = [-140,-129,55,60]
        flagregion.spatial_region = source_region
        flagregions.append(flagregion)

        # Westport Variable Region - 45sec to 9sec to 3sec:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Westport_45-9-3sec'
        flagregion.minlevel = 4
        flagregion.maxlevel = 6
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.29, -123.655, 46.325, 47.159]
        flagregions.append(flagregion)

        # Grays Harbor Region - allow  1" grids, require 9sec at least:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Grays_9-3-1sec'
        flagregion.minlevel = 6
        flagregion.maxlevel = 7
        flagregion.t1 = tstart_finestgrid
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.199, -123.809, 46.8, 47.145]
        flagregions.append(flagregion)

        # fixedgrid Region - require  3" grids, allow 1 and  1/3":
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'fixedgrid_3-1-13sec'
        flagregion.minlevel = 6
        flagregion.maxlevel = 8
        flagregion.t1 = tstart_finestgrid
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.185,-123.935,46.785,46.935]
        flagregions.append(flagregion)



    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    #rundata.gaugedata.gauges.append([1,-122.4, 47.781, 0., 1.e10])

    asce_gagues_file = '/Users/rjl/git/CopesHubTsunamis/info/asce_values.txt'
    asce_gauges = np.loadtxt(asce_gagues_file, skiprows=1)
    for k in range(0,len(asce_gauges),10):
        gaugeno = int(asce_gauges[k,0])
        gx = float(asce_gauges[k,1])
        gy = float(asce_gauges[k,2])
        if gy >= 40:
            rundata.gaugedata.gauges.append([gaugeno, gx, gy, 0, 1e9])
    VI_gagues_file = '/Users/rjl/git/CopesHubTsunamis/info/gaugesVI.txt'
    VI_gauges = np.loadtxt(VI_gagues_file, skiprows=1)
    for k in range(0,len(VI_gauges),1):
        gaugeno = int(VI_gauges[k,0])
        gx = float(VI_gauges[k,1])
        gy = float(VI_gauges[k,2])
        rundata.gaugedata.gauges.append([gaugeno, gx, gy, 0, 1e9])


    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    # Points on a uniform 2d grid:
    dx_fine = 0.0125  # grid resolution at finest level

    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid
    fg.x1 = -120. + dx_fine/2.  # specify pts to align with FV cell centers
    fg.x2 = -60. - dx_fine/2.
    fg.y1 = -60. + dx_fine/2.
    fg.y2 = 0. - dx_fine/2.
    fg.x1 = -130. + dx_fine/2.  # specify pts to align with FV cell centers
    fg.x2 = -122. - dx_fine/2.
    fg.y1 = 38.5 + dx_fine/2.
    fg.y2 = 50.5 - dx_fine/2.
    fg.dx = dx_fine
    fg.dy = dx_fine
    fg.tstart_max =  0.      # when to start monitoring max values
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 30.         # target time (sec) increment between updating
                              # max values
    fg.min_level_check = 4    # which levels to monitor max on
    fg.arrival_tol = 1.e-1    # tolerance for flagging arrival

    fg.interp_method = 0      # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)    # written to fgmax_grids.data


    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2       # will specify a 2d grid of points
    fgout.output_format = 'binary32'
    fgout.nx = 8*60
    fgout.ny = 12*60
    fgout.x1 = -130.  # specify edges (fgout pts will be cell centers)
    fgout.x2 = -122.
    fgout.y1 = 38.5
    fgout.y2 = 50.5
    fgout.tstart = 0.
    fgout.tend = 2*3600
    fgout.nout = int(np.floor(fgout.tend)/30) + 1
    fgout_grids.append(fgout)    # written to fgout_grids.data


    if 0:
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 2
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        fgout.nx = 990
        fgout.ny = 541
        fgout.x1 = -124.3  # specify edges (fgout pts will be cell centers)
        fgout.x2 = -123.75
        fgout.y1 = 46.8
        fgout.y2 = 47.1
        fgout.tstart = 0.
        fgout.tend = 60*60.
        fgout.nout = 121
        fgout_grids.append(fgout)    # written to fgout_grids.data



    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
    kmltools.make_input_data_kmls(rundata)
