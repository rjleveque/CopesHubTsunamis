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
region_tools  = fullpath_import(f'{common_code_dir}/region_tools.py')
kmltools  = fullpath_import(f'{common_code_dir}/kmltools.py')

HOME = os.environ['HOME']
# directory for shared topo and dtopo files:
CHTshare = CHT.replace(HOME,'/work2/04137/rjl/CHTshare')

# if some values in .data files show up as e.g. np.float64(3600.0)
# this will restore old behavior and just print 3600.0:
#np.set_printoptions(legacy="1.25")  # Loyce testing by commenting out
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
        #clawdata.tfinal = 10.     #test run for 10 seconds
        clawdata.tfinal = 10*3600. 
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
    try:
        geo_data.manning_coefficient = params.manning_coefficient
    except:
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
    #least likeable to more likeable
    rundata.topo_data.override_order = False

    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]

    # 15-second etopo22:
    topofiles.append([3, f'{topo_dir}/etopo22_15s_-137_-121_37_55.asc'])

    # 3 arcsec topo:
    topofiles.append([3, f'{topo_dir}/nw_pacific_3sec_WA_100m_gauges.asc'])

    # 2 arcsec topo:
    topofiles.append([3, f'{topo_dir}/astoriaS_2s_mhw.asc'])
    topofiles.append([3, f'{topo_dir}/GH_tiles_2021_2s.asc'])

    # 1/3 arcsec topo:
    topofiles.append([3, f'{topo_dir}/GH_tiles_2021_13s.asc'])

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
    flagregion.spatial_region = [-126.6,-123.6,46.15,47.6185]
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
    gauges_region = [-125.9,-123.75,46.6,47.227]
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
    flagregion.spatial_region = [-124.295,-123.75,46.6,47.227]
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
    flagregion.spatial_region = [-124.245,-123.95, 46.9,47.145]
    flagregions.append(flagregion)
    
    # North Shore of Grays Harbor, Humptulips Valley
    # Level 8 is 1/3 sec
    # Note that this region contains the fgmax_extent, but is larger
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_13sec_NS1'
    flagregion.minlevel = 8
    flagregion.maxlevel = 8
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath('RuledRectangle_NS1.data')
    flagregions.append(flagregion)

    slu_NS1=np.array([[-124.074,47.02,47.107],[-124.0412,47.02,47.107],
                      [-124.016,47.02,47.07],[-123.995,47.02,47.07]])
    rr_NS1 = region_tools.RuledRectangle(slu=slu_NS1)
    rr_NS1.ixy = 'x'     #s refers to x, the longitude, lower and upper in lat y
    rr_NS1.method = 1  #connect linearly
    rr_NS1_name = 'RuledRectangle_NS1'
    rr_NS1.write(rr_NS1_name + '.data')

    # G80 includes gauge 80
    # Level 8 is 1/3 sec
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_13sec_G80'
    flagregion.minlevel = 8
    flagregion.maxlevel = 8
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.03,-123.995,46.992,47.02]
    flagregions.append(flagregion)

    # North Shore of Grays Harbor, Ocean past Hogan's Corner and road out
    # Level 8 is 1/3 sec
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_13sec_NS2'
    flagregion.minlevel = 8
    flagregion.maxlevel = 8
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath('RuledRectangle_NS2.data')
    flagregions.append(flagregion)

    slu_NS2=np.array([[ -124.185,47.02,47.053],[-124.1625,47.02,47.053],
                      [-124.1225,47.02,47.06],
                      [-124.074,47.02,47.06]])
    rr_NS2 = region_tools.RuledRectangle(slu=slu_NS2)
    rr_NS2.ixy = 'x'     #s refers to x, the longitude, lower and upper in lat y
    rr_NS2.method = 1  #connect linearly
    rr_NS2_name = 'RuledRectangle_NS2'
    rr_NS2.write(rr_NS2_name + '.data')

    # ---------------
    # Gauges:
    # ---------------
    #For Bridges
    gauges = []
    if 1: #from NS1 region
        gauges.append([61, -124.04981481, 47.05657407, 0, 1e+09])
        gauges.append([74, -124.04018519, 47.04805556, 0, 1e+09])
        gauges.append([75, -124.04287037, 47.05111111, 0, 1e+09])
        gauges.append([81, -124.02546296, 47.03111111, 0, 1e+09])
        gauges.append([221, -124.03851852, 47.07268519, 0, 1e+09])
        gauges.append([222, -124.04277778, 47.07203704, 0, 1e+09])
        gauges.append([223, -124.04361111, 47.07194444, 0, 1e+09])
        gauges.append([229, -124.05083333, 47.07314815, 0, 1e+09])
        gauges.append([262, -124.04027778, 47.06462963, 0, 1e+09])
        gauges.append([292, -124.05703704, 47.04537037, 0, 1e+09])
        gauges.append([293, -124.05898148, 47.04472222, 0, 1e+09])
        gauges.append([331, -124.04472222, 47.10129630, 0, 1e+09])
        gauges.append([336, -124.06287037, 47.10296296, 0, 1e+09])
        gauges.append([355, -124.05240741, 47.07453704, 0, 1e+09])
        gauges.append([356, -124.05935185, 47.09740741, 0, 1e+09])
        gauges.append([434, -124.05018519, 47.07194444, 0, 1e+09])
        gauges.append([438, -124.00379630, 47.05944444, 0, 1e+09])
    if 1:  #from G80 region
        gauges.append([80, -124.00175926, 47.00416667, 0, 1e+09])
    rundata.gaugedata.gauges = gauges

    # ---------------
    # fgout grids:
    # ---------------

    # ---------------
    # fgmax grids:
    # ---------------


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

    # run this as script via 'python setrun_case.py`
    # to make data without any dtopofile, to check other inputs:
    rundata = setrun('geoclaw', case={})
    rundata.write()

    # To create kml files of inputs:
    kmltools.make_input_data_kmls(rundata)
