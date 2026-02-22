"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

Note: this version has case as a parameter to setrun, a dictionary used to pass
in values that vary from case to case for doing a parameter sweep.
"""
import sys,os,datetime
import numpy as np

from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgmax_tools, fgout_tools
from clawpack.geoclaw.data import ForceDry
from clawpack.clawutil.util import fullpath_import
import pathlib

# top level directory for this project:
try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

common_code_dir = f'{CHT}/common_code'
restart_tools = fullpath_import(f'{common_code_dir}/restart_tools.py')

HOME = os.environ['HOME']
# directory for shared topo and dtopo files:
CHTshare = CHT.replace(HOME,'/work2/04137/rjl/CHTshare')

# if some values in .data files show up as e.g. np.float64(3600.0)
# this will restore old behavior and just print 3600.0:
# np.set_printoptions(legacy="1.25")
# fixed in master after v5.12.0


print('\n========================== setrun.py ==========================')
print('Start date & time: ', datetime.datetime.now())


use_bouss_version = False  # True ==> requires code compiled with PETSc etc.
                           # might still solve SWE if bouss_equations = 0 below

topo_dir = f'{CHT}/topo/topofiles'

# don't need dtopo_dir here since dtopofile path will be passed in
# to setrun in case dictionary by runclaw_makeplots_dtopos.py
#dtopo_dir = f'{CHT}/dtopo/CSZ_groundmotions/dtopo30sec/dtopofiles'

# location for big files for different computer environments:
# output and plots will be sent to scratch_dir

this_dir = os.getcwd()

if 'rjl/git' in this_dir:
    computer = 'rjl-laptop'
    scratch_dir = this_dir.replace('rjl/git/CopesHubTsunamis/geoclaw_runs', \
                                   'rjl/scratch/CHT_runs')
    topo_dir = '/Users/rjl/topo/topofiles'

elif '/mmfs1/home' in this_dir:
    computer = 'hyak'
    scratch_dir = this_dir.replace('/mmfs1/home', '/gscratch/tsunami')
    topo_dir = topo_dir.replace('/mmfs1/home', '/gscratch/tsunami')

elif '/home1' in this_dir:
    computer = 'tacc'
    scratch_dir = this_dir.replace('/home1', '/scratch')
    #topo_dir = topo_dir.replace('/home1', '/scratch')
    topo_dir = f'{CHTshare}/topo/topofiles'

else:
    computer = 'unknown'
    scratch_dir = this_dir



print('topo_dir is:  ',topo_dir)

RRdir = f'{CHT}/topo/regions'  # for flagregion Ruled Rectangles



#------------------------------
def setrun(claw_pkg='geoclaw', case={}):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    ----------------
    case parameters:
    ----------------

    This version has a parameter `case` that is expected to be a dictionary
    with the key 'dtopofile' giving the full path to the dtopofile for this
    case. (If not, it's assumed that no dtopofile should be used.)

    dtopo_type is assumed to be 3, unless case['dtopo_type'] is set otherwise

    case['outdir'] is generally set but is only needed for a restart, in
    which case it is used to determine which checkpoint file to restart from.

    """

    from clawpack.clawutil import data
    from pathlib import Path

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)



    # full path to dtopofile to use for this case:
    dtopofile = case.get('dtopofile', None) # returns None if not set

    dtopo_type = case.get('dtopo_type', 3)  # assume dtopo_type=3 if not set

    if dtopofile is not None:
        event = Path(dtopofile).stem  # drop path and .dtt3 extension
    else:
        event = 'NO_DTOPO'

    # We will set the initial time step dt_initial and
    # dtopo_data.dt_max_dtopo differently for instant vs kinematic ruptures:

    instant = 'instant' in event  # part of filename for static ruptures
    if instant:
        # instantaneous static rupture
        dt_max_dtopo = 0.25
        print(f'Assuming event {event} is instantaneous displacement')
    else:
        # kinematic rupture
        dt_max_dtopo = 5.0
        print(f'Assuming event {event} is kinematic displacement')
    print(f'    dt_max_dtopo = {dt_max_dtopo:.1f} seconds')

    restart = False
    restart_file = ''

    if 0:
        # restart from previous run, if checkpt file is available:
        # (note: requires case['outdir'] to be set by calling program)
        restart_file = restart_tools.find_last_checkpt(case['outdir'])

    if restart_file != '':
        restart = True
        restart_time = restart_tools.time(restart_file)
        print(f'Will restart from time t = {restart_time}')


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

    sec16 = 1./(6*3600.)  # one sixth arcsecond
    clawdata.lower[0] = -136.75 - sec16       # west longitude
    clawdata.upper[0] = -122.75 - sec16       # east longitude

    clawdata.lower[1] = 38.0 - sec16         # south latitude
    clawdata.upper[1] = 54.0 - sec16         # north latitude

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] =  14*15
    clawdata.num_cells[1] =  16*15

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    if use_bouss_version:
        clawdata.num_eqn = 5
    else:
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

    clawdata.restart = False
    clawdata.restart_file = 'fort.chkbbbbb'

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
        # Here we run for 120 seconds and produce NO time frame output:
        clawdata.num_output_times = 0
        #clawdata.tfinal  = 20*60    #20 minute test on laptop
        #clawdata.tfinal = 10*3600. 
        clawdata.tfinal = 2*3600. 
        ##################clawdata.tfinal = 10*3600.
        clawdata.output_t0 = False

    elif clawdata.output_style == 2:
        # Specify a list of output times. Params will have the times.
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
    clawdata.dt_initial = dt_max_dtopo

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.85
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
        clawdata.checkpt_times = 5*3600*np.arange(1,4,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 100


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata
    amrdata.max1d = 200
    amrdata.memsize = 2**27 - 1

    # max number of refinement levels:
    #amrdata.amr_levels_max = 6    # for 1/3" resolution
    #amrdata.amr_levels_max = 4    # for 3"
    #amrdata.amr_levels_max = 1     # for short test

    amrdata.amr_levels_max = 8    # for real test

    # List of refinement ratios at each level (length at least mxnest-1)

    # Set up for 9 levels below, using only 8 for now.
    # dx = dy = 4', 2', 24", 12", 6", 3", 1", 1/3", 1/9"
    refinement_ratios = [2,5,2,2,2,3,3,3]
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
    geo_data.speed_limit = 20.
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.05

    # ---------------
    # Topo:
    # ---------------
    # == settopo.data values ==

    #Use the topo files prioritized as listed below from
    #least likeable to more likeable -- dont use has problems
    #Always set to False below
    rundata.topo_data.override_order = False 

    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]

    # 15-second etopo22:
    topofiles.append([3, f'{topo_dir}/etopo22_15s_-137_-121_37_55.asc'])

    # 3 arcsec topo:
    topofiles.append([3, f'{topo_dir}/nw_pacific_3sec_WA_100m_gauges.asc'])

    # 1 arcsec topo:
    topofiles.append([3, f'{topo_dir}/n46x75_w124x00_1s.asc'])
    topofiles.append([3, f'{topo_dir}/n46x75_w124x25_1s.asc'])
    topofiles.append([3, f'{topo_dir}/n47x00_w124x00_1s.asc'])
    topofiles.append([3, f'{topo_dir}/n47x00_w124x25_1s.asc'])
    
    # 1/3 arcsec topo:
    topofiles.append([3, f'{topo_dir}/WillapaN_13s.asc'])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    if dtopofile is not None:
        dtopo_data.dtopofiles = [[dtopo_type, dtopofile]] # from case dictionary

    # maximum time step to use on level 1 while rupturing:
    dtopo_data.dt_max_dtopo = dt_max_dtopo


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

    # Computational domain Variable Region - 4 minutes to 24":
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [clawdata.lower[0]-0.2,
                                 clawdata.upper[0]+0.2,
                                 clawdata.lower[1]-0.2,
                                 clawdata.upper[1]+0.2]
    flagregions.append(flagregion)

    # Region12sec - fixed at 12 sec :
    # Level 4 is 12 sec
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_12sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-126.6,-123.6,45.92,47.39]
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 24", 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51b_24sec_12sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51b.data')
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 24", 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46b_24sec_12sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46b.data')
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_12sec_6sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 5
    flagregion.t1 = 0.0
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle for now
    gauges_region = [-125.9,-123.75,46.27,47.08]
    flagregion.spatial_region = gauges_region
    flagregions.append(flagregion)

    # 6sec to 3sec region
    # Level 6 is 3 sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_6sec_3sec'
    flagregion.minlevel = 5
    flagregion.maxlevel = 6
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.295,-123.75,46.27,47.08]
    flagregions.append(flagregion)

    # 3sec to 1sec region
    # Level 7 is 1 sec
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_3sec_1sec'
    flagregion.minlevel = 6
    flagregion.maxlevel = 7
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.175,-123.87, 46.65,46.755]
    flagregions.append(flagregion)

    # SBIT 1/3sec region 
    # Level 8 is 1/3 sec
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_13sec_SBIT'
    flagregion.minlevel = 8
    flagregion.maxlevel = 8
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.1,-123.965,46.680,46.735]
    flagregions.append(flagregion)

    # ---------------
    # Gauges:
    # ---------------
    #
    #For Gauges for SBIT, put here later
    gauges = []
    if 0:  #Desired at the moment, might need to center them
        #Eagle Hill Rd high ground entrance
        gauges.append([1, -124.029176, 46.728031, 0, 1e+09])

        #Eagle Hill Rd and SR 105
        gauges.append([2, -124.027977, 46.726885, 0, 1e+09])

        #Tokeland Rd and SR 105
        gauges.append([3, -124.020344, 46.72461, 0, 1e+09])

        #Shoalwater Bay Gynasium
        gauges.append([4, -124.013834, 46.721922, 0, 1e+09])

        #Shoalwater Bay Gynasium Trail Middle
        gauges.append([5, -124.014507, 46.723688, 0, 1e+09])

        #Shoalwater Bay Gynasium Trail and SR 105
        gauges.append([6, -124.013470, 46.725850, 0, 1e+09])

        #Trail to Teal Duck Slough, closer to SR 105
        gauges.append([7, -124.010387, 46.727853, 0, 1e+09])

        #Trail to High Ground past Duck Slough Road and SR 105
        gauges.append([8, -124.004113, 46.727952, 0, 1e+09])

        #High Ground past Duck Slough Road (clearcut)
        gauges.append([9, -124.005254, 46.728384, 0, 1e+09])

        #Shoalwater Bay Tribal Court
        gauges.append([10, -124.016535, 46.721395, 0, 1e+09])

        #Shoalwater Bay Food Bank
        gauges.append([11, -124.0152, 46.721125, 0, 1e+09])

        #Tradewinds on the Bay
        gauges.append([12, -124.01304, 46.717688, 0, 1e+09])

        #Auntie Lee and Blackberry Lane entrance
        gauges.append([13, -123.990825, 46.710251, 0, 1e+09])

        #South Beach Regional fire Authority
        gauges.append([14, -123.995923, 46.711095, 0, 1e+09])

        #Shoalwater Museum
        gauges.append([15, -124.021153, 46.724253, 0, 1e+09])


    if 1:
        #Centered version of the desired 14 gauges above
        gauges.append([1, -124.02916667, 46.72805556, 0, 1e+09])
        gauges.append([2, -124.02796296, 46.72685185, 0, 1e+09])
        gauges.append([3, -124.02037037, 46.72462963, 0, 1e+09])
        gauges.append([4, -124.01379630, 46.72194444, 0, 1e+09])
        gauges.append([5, -124.01453704, 46.72370370, 0, 1e+09])
        gauges.append([6, -124.01342593, 46.72583333, 0, 1e+09])
        gauges.append([7, -124.01037037, 46.72787037, 0, 1e+09])
        gauges.append([8, -124.00407407, 46.72796296, 0, 1e+09])
        gauges.append([9, -124.00527778, 46.72842593, 0, 1e+09])
        gauges.append([10, -124.01657407, 46.72138889, 0, 1e+09])
        gauges.append([11, -124.01518519, 46.72111111, 0, 1e+09])
        gauges.append([12, -124.01305556, 46.71768519, 0, 1e+09])
        gauges.append([13, -123.99083333, 46.71027778, 0, 1e+09])
        gauges.append([14, -123.99592593, 46.71111111, 0, 1e+09])
        gauges.append([15, -124.02111111, 46.72425926, 0, 1e+09])

    #For Bridges
    if 0:   #Aberdeen Bridges
        gauges.append([11, -123.81157407, 46.97601852, 0, 1e+09])
        gauges.append([33, -123.79111111, 46.97833333, 0, 1e+09])
        gauges.append([43, -123.81185185, 46.97712963, 0, 1e+09])
        gauges.append([56, -123.77611111, 46.96000000, 0, 1e+09])
        gauges.append([68, -123.80861111, 46.97212963, 0, 1e+09])
        gauges.append([137, -123.81212963, 46.97407407, 0, 1e+09])
        gauges.append([138, -123.81296296, 46.97342593, 0, 1e+09])
        gauges.append([354, -123.80074074, 46.99842593, 0, 1e+09])
        gauges.append([367, -123.78148148, 46.97962963, 0, 1e+09])
        gauges.append([368, -123.80518519, 46.98462963, 0, 1e+09])
        gauges.append([390, -123.80129630, 46.97694444, 0, 1e+09])
        gauges.append([509, -123.77787037, 46.95972222, 0, 1e+09])
        gauges.append([558, -123.77037037, 46.97490741, 0, 1e+09])
        gauges.append([561, -123.81009259, 46.97481481, 0, 1e+09])
        gauges.append([562, -123.81009259, 46.97481481, 0, 1e+09])
    rundata.gaugedata.gauges = gauges

    # ---------------
    # fgmax grids:
    # ---------------
    # NEW STYLE STARTING IN v5.7.0
    # ------------------------------
    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 5

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    #Set fgmax_extent: want the boundaries to be cell centers, so
    #can be an even integer plus some multiple of 1/3sec.  This will be true
    #if the decimal part is a multiple of .1, or .01, .005, or .0025 for
    #example based on how the computational domain was shifted.  This

    fgmax_extent=[-124.0375,-123.99,46.7075,46.735]
    if 1:
        # Points on a uniform 2d grid:
        dx_fine = 1./(3*3600.)  # grid resolution at finest level
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = fgmax_extent[0] #+ dx_fine/2.
        fg.x2 = fgmax_extent[1] #- dx_fine/2.
        fg.y1 = fgmax_extent[2] #+ dx_fine/2.
        fg.y2 = fgmax_extent[3] #- dx_fine/2.
        fg.dx = dx_fine
        fg.tstart_max =  clawdata.t0 + 120   # when to start monitoring max values
        fg.tend_max = 1.e10         # when to stop monitoring max values
        fg.dt_check = 5.           # target time (sec) increment between updating
                                    # max values
                            # which levels to monitor max on
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-2      # tolerance for flagging arrival

        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data

    # ---------------
    # fgout grids:
    # ---------------

    if 1:
        ###  HERE IS THE FGOUT STUFF TO EDIT (if plot_eta is set to True in
        ###  compare_gauge_max_withasce.py postprocessing, it uses fgout)
        # == fgout_grids.data values ==
        # NEW IN v5.9.0
        # Set rundata.fgout_data.fgout_grids to be a list of
        # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:

        #This is the region where 1/3 computation was
        #flagregion.spatial_region = [-124.1,-123.965,46.680,46.735]
        one_sixth = 1.0/(6.0*3600.)

        fgout_grids = rundata.fgout_data.fgout_grids  #empty list initially
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 1
        fgout.point_style = 2  #a 2d grid of points
        fgout.output_format = 'binary32'
        fgout.nx = 1459
        fgout.ny = 595
        fgout.x1 = -124.1 - one_sixth
        fgout.x2 = -123.965 + one_sixth
        fgout.y1 = 46.68 - one_sixth
        fgout.y2 = 46.735 + one_sixth
        fgout.tstart = 0.0
        fgout.tend = 2*3600
        fgout.nout = int(np.round(fgout.tend/20)) + 1
        fgout_grids.append(fgout)    # written to fgout_grids.data

    if use_bouss_version:
        # To use Boussinesq solver, add bouss_data parameters here
        # Also make sure to use the correct Makefile pointing to bouss version
        # and set clawdata.num_eqn = 5

        from clawpack.geoclaw.data import BoussData
        rundata.add_data(BoussData(),'bouss_data')

        rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN
        rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss
        rundata.bouss_data.bouss_max_level = 7    # finest level to apply bouss
        rundata.bouss_data.bouss_min_depth = 10.  # depth to switch to SWE
        rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc
        rundata.bouss_data.bouss_tstart = 60.      # time to switch from SWE


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

    # run this as script via 'python setrun_case.py`
    # to make data without any dtopofile, to check other inputs:
    rundata = setrun('geoclaw', case={})
    rundata.write()

    # To create kml files of inputs:
    kmltools.make_input_data_kmls(rundata)
