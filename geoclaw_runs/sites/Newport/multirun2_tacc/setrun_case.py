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


# if some values in .data files show up as e.g. np.float64(3600.0)
# this will restore old behavior and just print 3600.0:
np.set_printoptions(legacy="1.25")
# fixed in master after v5.12.0


print('\n========================== setrun.py ==========================')
print('Start date & time: ', datetime.datetime.now())


use_bouss_version = False  # True ==> requires code compiled with PETSc etc.
                           # might still solve SWE if bouss_equations = 0 below


# Set these directories where input data is found:

rundir = os.getcwd()
print('rundir = %s' % rundir)

topo_dir = f'{CHT}/topo/topofiles'

# for hyak cluster:
#topo_dir = topo_dir.replace('/mmfs1/home', '/gscratch/tsunami')

# for TACC:
topo_dir = topo_dir.replace('/home1', '/scratch')

# laptop:
#topo_dir = '/Users/rjl/topo/topofiles'

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

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # The values below are expected to be in case dictionary,
    # and may vary from case to case:

    dtopofiles = case['dtopofiles']

    if len(dtopofiles) > 0:
        dtopofile = pathlib.Path(dtopofiles[0][1])
        event = dtopofile.stem
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

    # Lower and upper edge of computational domain:
    #clawdata.lower[0] = -128.      # west longitude
    #clawdata.lower[1] = 38.5       # south latitude
    #clawdata.upper[1] = 53.5       # north latitude
    #clawdata.lower[1] = 44.5       # south latitude
    #clawdata.upper[1] = 47.5       # north latitude
    #clawdata.num_cells[0] = 6*15
    #clawdata.num_cells[1] = 3*15

    one_sixth = 1.0/(6.0*3600.)
    one_sixth = 0.  # FIX
    clawdata.lower[0] = -135. - one_sixth      # west longitude
    clawdata.upper[0] = -122. - one_sixth      # east longitude
    clawdata.lower[1] = 38.5  - one_sixth      # south latitude
    clawdata.upper[1] = 51.5  - one_sixth      # north latitude

    clawdata.num_cells[0] = 13*15
    clawdata.num_cells[1] = 13*15


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

        clawdata.num_output_times = 0      # no frame output
        clawdata.tfinal = 2.0*3600.
        clawdata.output_t0 = False        # output at initial (or restart) time?

        if 0:
            clawdata.num_output_times = 1
            clawdata.tfinal = 30.              # SHORT TEST
            clawdata.output_t0 = True          # output at initial (or restart) time?

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
    clawdata.dt_initial = dt_max_dtopo

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.8
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
        clawdata.checkpt_times = 1800.*np.arange(1,11,1)

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
    amrdata.amr_levels_max = 6   # for 1/3" resolution
    #amrdata.amr_levels_max = 4    # for 3"

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 1deg, 6', 3', 45", 9", 3", 1", 1/3"":
    # refinement_ratios = [10,2,4,5,3,3,3]

    # dx = dy = 1deg, 6', 1', 15"
    #refinement_ratios = [10,6,4]

    # dx = dy = 1deg, 6', 30", 7.5"
    #refinement_ratios = [10,12,4]

    # dx = dy = 1deg, 6', 30", 15", 5"
    #refinement_ratios = [10,12,2,3]

    # dx = dy = 4', 1', 30", 15", 5", 1", 1/3", 1/9"
    #refinement_ratios = [4,2,2,3,5,3,3]

    # dx = dy = 4', 1', 30", 10", 5", 1", 1/3", 1/9"
    #refinement_ratios = [4,2,3,2,5,3,3]

    # dx = dy = 4', 2', 24", 12", 6", 3", 1", 1/3", 1/9"
    #refinement_ratios = [2,5,2,2,2,3,3,3]

    # dx = dy = 4', 1', 12", 3", 1", 1/3", 1/9"
    refinement_ratios = [4,5,4,3,3,3]


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
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]

    # 15-second etopo22:
    topofiles.append([3, f'{topo_dir}/etopo22_15s_-137_-121_37_55.asc'])

    #3-second topo:
    topofiles.append([3, f'{topo_dir}/crm_vol8_3sec_cropped_44-47.asc'])

    # 1/3-second topo:
    topofiles.append([3, f'{topo_dir}/NewportE_13s_mhw.asc'])
    topofiles.append([3, f'{topo_dir}/Newport_13s_mhw.asc'])

    # 1/9-second topo:
    topofiles.append([3, f'{topo_dir}/YaquinaBay_19s_mhw.asc'])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    dtopo_data.dtopofiles = dtopofiles # from case dictionary

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

    # here we just use flagregions:
    flagregions = rundata.flagregiondata.flagregions  # initialized to []

    # Computational domain Variable Region
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [clawdata.lower[0]-0.1,
                                 clawdata.upper[0]+0.1,
                                 clawdata.lower[1]-0.1,
                                 clawdata.upper[1]+0.1]
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51b_12sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    #flagregion.t2 = tmax_dtopo_region
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51b.data')
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46b_12sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 3
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46b.data')
    flagregions.append(flagregion)



    if 0:
        # No longer have 6" level

        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_6sec'
        flagregion.minlevel = 5
        flagregion.maxlevel = 5
        flagregion.t1 = 0.0
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle for now

        ## for Ocean Shores, encompasses gauges 94 to 137
        #latitudes below are for 93 in north to 138 in south
        #gauges_region = [-125.9,-124.1,46.66,47.287]

        ## for Westport, encompasses gauges 98 to 143
        #latitudes below are for 97 in north to 144 in south
        #gauges_region = [-125.9,-124.1,46.587,47.227]

        ## for Seaside, encompasses gauges 166 to 208
        #gauges_region = [-125.65,-123.76,45.69,46.43]
        #gauges_region = [-126.00,-123.76,45.69,46.43]

        ## for Newport, encompasses  +-32.2km of the gauges on land
        gauges_region = [-125.75,-123.9125,44.29,44.96]

        flagregion.spatial_region = gauges_region
        flagregions.append(flagregion)


    if 1: #For Newport

        # Old Region 2" - fixed at 3" sec. :
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_3sec'
        flagregion.minlevel = 3
        flagregion.maxlevel = 4
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.25,-123.9150,44.535,44.715]
        flagregions.append(flagregion)

        # Region 1" - fixed at 1" sec
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_1sec'
        flagregion.minlevel = 3
        flagregion.maxlevel = 5
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.1495,-123.916,44.55,44.6795]
        flagregions.append(flagregion)

        # Region 1/3" - fixed at 1/3"
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_onethird'
        flagregion.minlevel = 3
        flagregion.maxlevel = 6
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.102,-123.9175,44.564,44.675]
        flagregions.append(flagregion)

    if 0:  #set up, but not using this flagregion at the moment
        # Level 7 is 1/9" sec, might need around marina?
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_oneninth'
        flagregion.minlevel = 3
        flagregion.maxlevel = 7
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.0895,-124.002,44.601,44.634]
        flagregions.append(flagregion)

    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    rundata.gaugedata.file_format = 'binary32' # save as single precision
    rundata.gaugedata.min_time_increment = 5 # seconds min between outputs

    if 1:
        # load list of virtual gauges to use, with columns
        #      gaugeno, x, y
        # x,y should be in decimal form, preferably cell centered on finest grid
        gauges_file = f'{CHT}/gauges/oregon_gauges_2024/VGListNewport.csv'
        gauges = np.loadtxt(gauges_file, delimiter=',')

        for k in range(gauges.shape[0]):
            gaugeno = int(gauges[k,0])
            xg = gauges[k,1]
            yg = gauges[k,2]
            rundata.gaugedata.gauges.append([gaugeno,xg,yg,0,1e9])

    # fgmax grids:

    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed

    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    if 1:
        # fgmax region covering Yaquina Bay and Newport, west of -124.
        fg = fgmax_tools.FGmaxGrid()
        fg.xy_fname = f'{CHT}/topo/fgmax_grids/fgmax_pts_Newport.data'
        fg.point_style = 4
        fg.tstart_max =  30.*60.  # when to start monitoring max values, at 20minutes
        fg.tend_max = 1.e10       # when to stop monitoring max values
        fg.dt_check = 5.          # target time (sec) increment between updating.
        # monitor on finest level:
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-1    # tolerance for flagging arrival
        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data

    if 1:
        # fgmax region covering Yaquina River, east of -124.
        fg = fgmax_tools.FGmaxGrid()
        fg.xy_fname = f'{CHT}/topo/fgmax_grids/fgmax_pts_NewportE.data'
        fg.point_style = 4
        fg.tstart_max =  30.*60.  # when to start monitoring max values, at 20minutes
        fg.tend_max = 1.e10       # when to stop monitoring max values
        fg.dt_check = 5.          # target time (sec) increment between updating.
        # monitor on finest level:
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-1    # tolerance for flagging arrival
        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data



    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    if use_bouss_version:
        q_out_vars = [1,2,3,6] # h,hu,hv,eta for bouss code
    else:
        q_out_vars = [1,2,3,4] # h,hu,hv,eta for shallow code

    if 0:
        # full domain (smaller for inundation run than for offshore_gauges)
        dx_fgout = 120./3600.  # degrees (two minute)
        dy_fgout = 120./3600.  # degrees
        dt_fgout = 30  # seconds
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 1
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        fgout.x1 = clawdata.lower[0]    # specify edges
        fgout.x2 = clawdata.upper[0]    # (fgout pts will be cell centers)
        fgout.y1 = clawdata.lower[1]    # edge of a cell
        fgout.y2 = clawdata.upper[1]    # edge of a cell
        fgout.nx = int(round((fgout.x2 - fgout.x1)/dx_fgout))
        fgout.ny = int(round((fgout.y2 - fgout.y1)/dy_fgout))
        fgout.tstart = 0.
        fgout.tend = clawdata.tfinal
        fgout.nout = int(np.floor((fgout.tend - fgout.tstart))/dt_fgout) + 1
        fgout.q_out_vars = q_out_vars
        fgout_grids.append(fgout)    # written to fgout_grids.data


    if 0:
        # small region for inset plot
        dx_fgout = 6/3600.
        dy_fgout = 12/3600.  # degrees
        dt_fgout = 30  # seconds
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 2
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        fgout.x1 = -124.145 - one_sixth  # specify edges (fgout pts will be cell centers)
        fgout.x2 = -124.055 - one_sixth
        fgout.y1 = 44.56 - one_sixth
        fgout.y2 = 44.67 - one_sixth
        fgout.nx = int(round((fgout.x2 - fgout.x1)/dx_fgout))
        fgout.ny = int(round((fgout.y2 - fgout.y1)/dy_fgout))
        fgout.tstart = 0.
        fgout.tend = clawdata.tfinal
        fgout.nout = int(np.floor((fgout.tend - fgout.tstart))/dt_fgout) + 1
        fgout.q_out_vars = q_out_vars
        fgout_grids.append(fgout)    # written to fgout_grids.data

    if 0:
        #Keep this inside the 1/3 computational grid which was
        #[-124.102,-123.9175,44.564,44.675]

        west13 = -124.1 - one_sixth; east13 = -123.9175 - one_sixth;
        south13 = 44.565 - one_sixth; north13 = 44.675 - one_sixth;
        fgout_extent = [west13,east13,south13,north13]
        dx_fgout = 1/3 * 1/3600.  # degrees
        dy_fgout = dx_fgout
        dt_fgout = 15  # seconds

        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 3
        fgout.point_style = 2  # uniform rectangular x-y grid
        fgout.output_format = 'binary32'
        fgout.x1 = fgout_extent[0]
        fgout.x2 = fgout_extent[1]
        fgout.y1 = fgout_extent[2]
        fgout.y2 = fgout_extent[3]
        fgout.nx = int(round((fgout.x2 - fgout.x1)/dx_fgout))
        fgout.ny = int(round((fgout.y2 - fgout.y1)/dy_fgout))
        fgout.tstart = 20*60.
        #fgout.tstart = 0.  # SHORT TEST
        fgout.tend = clawdata.tfinal
        fgout.nout = int(round(((fgout.tend - fgout.tstart)/dt_fgout))) + 1
        fgout.nout = max(fgout.nout, 0)
        fgout.q_out_vars = q_out_vars
        fgout_grids.append(fgout)


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
    rundata = setrun('geoclaw', case={'dtopofiles':[]})
    rundata.write()

    # To create kml files of inputs:
    kmltools.make_input_data_kmls(rundata)
