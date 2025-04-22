"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""
import sys
import os
import numpy as np
import datetime
#print('+++ cwd: ',os.getcwd())
sys.path.insert(0,'.')
#print('+++ sys.path: ',sys.path)


from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgmax_tools, fgout_tools
from clawpack.geoclaw.data import ForceDry

# top level directory for this project:
root_dir     = os.environ['CHT']
rundir = os.getcwd()

topo_dir = root_dir + '/topo/topofiles'
dtopo_dir = root_dir + '/dtopo/CSZ_groundmotions/dtopofiles'
RRdir = root_dir + '/topo/regions'

# for hyak cluster:
topo_dir = topo_dir.replace('/mmfs1/home', '/gscratch/tsunami')
dtopo_dir = dtopo_dir.replace('/mmfs1/home', '/gscratch/tsunami')

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

print('\n==================================== setrun.py ====================================')
print('Start date & time: ', datetime.datetime.now())

# Set these directories where input data is found:

rundir = os.getcwd()
print('rundir = %s' % rundir)


dtopofile = dtopo_dir + '/buried_locking_str10-shallow'
assert os.path.isfile(dtopofile), '*** did not find dtopofile = %s' \
        % dtopofile

input_dir            = root_dir + '/topo/input_files'
gauges_dir           = root_dir + '/info/gauges'

print('topo_dir is:  ',topo_dir)
print('input_dir is:  ',input_dir)
print('dtopofile is:  ',dtopofile)
print('gauges_dir is:  ',gauges_dir)


if 0:

    gauges_info_csv = root_dir + \
        '/gauges/GraysHarborPacificCountyBridges/GraysHarborPacificCounty.csv'
    print(' ')
    print('gauges_info_csv file was: ',gauges_info_csv)
    print(' ')
    gauges_info = open(gauges_info_csv).readlines()

    # extent of gauges to consider in this run:
    yg1 = 46.77
    yg2 = 47.36

    gauge_x = {}
    gauge_y = {}
    gauge_flag = {}
    for line in gauges_info[1:]:
        line = line.split(',')
        gauge_no = int(line[0])
        gauge_x[gauge_no] = float(line[2])
        gauge_y[gauge_no] = float(line[1])
        gauge_flag[gauge_no] = int(line[5])
        if (gauge_y[gauge_no] < yg1) or (gauge_y[gauge_no] > yg2):
            # ignore gauges outside this extent:
            gauge_flag[gauge_no] = 0




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

    # shift so that cell centers on finest grid align with DEM:
    sec16 = 1./(6*3600.)  # one sixth arcsecond

    # Use (-127.6,-123.0,39.6,49.6) shifted by sec16
    #BIG DOMAIN uses -137,-123,38,54
    #For now, don't have topo out to -137 -sec16, so change a bit

    #clawdata.lower[0] = -137 - sec16            # west longitude
    clawdata.lower[0] = -136.75 - sec16          # west longitude
    #clawdata.upper[0] = -123 - sec16            # east longitude
    clawdata.upper[0] = -122.75 - sec16          # east longitude

    clawdata.lower[1] = 38 - sec16         # south latitude
    clawdata.upper[1] = 54 - sec16         # north latitude

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] =  14*15  #210
    clawdata.num_cells[1] =  16*15  #240

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

    t0 = 0            # Start at time 0
    clawdata.t0 = t0  # Start time in seconds

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False     # True to restart from prior results
    clawdata.restart_file = ''   # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 2

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 20
        clawdata.tfinal = 5*3600.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        tfinal = 1.5*3600.
        dtout = 5*60.
        clawdata.output_times = [0.,60.,5*60.,10*60.] + \
                                list(np.arange(30*60., tfinal+1, dtout))

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True


    clawdata.output_format = 'binary'        # 'ascii' or 'binary'
    clawdata.output_q_components = 'all'     # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False     # output aux arrays each frame

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

    # Desired Courant number if variable dt used
    clawdata.cfl_desired = 0.74
    # max Courant number to allow without retaking step with a smaller dt:
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
    # clawdata.source_split = 'godunov'
    clawdata.source_split = 1

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
        clawdata.checkpt_times = 3600.*np.arange(1,6,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 8

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 4', 1', 30", 6", 2", 1", 1/3"":
    #amrdata.refinement_ratios_x = [4,2,5,3,2,3]
    #amrdata.refinement_ratios_y = [4,2,5,3,2,3]
    #amrdata.refinement_ratios_t = [4,2,5,3,2,3]

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
    amrdata.clustering_cutoff = 0.7

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
    geo_data.friction_depth = 200

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.2

    # ---------------
    # TOPO:
    # ---------------
    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]

    topofiles = topo_data.topofiles

    # 15-sec topo:
    topo_file = os.path.join(topo_dir, 'etopo22_15s_-137_-121_37_55.asc')
    topofiles.append([3, topo_file])

    if 1:
        # 3 arcsec topo
        #topo_file = os.path.join(topo_dir, 'nw_pacific_3sec_SW_WA.asc')
        #topo_file = os.path.join(topo_dir, 'nw_pacific_3sec_WA_100m_gauges.asc')
        topo_file = os.path.join(topo_dir, 'crm_vol8_3sec_cropped_44-47.asc')

        topofiles.append([3, topo_file])

    if 1:
        # 2 arc sec topo:
        topo_file = os.path.join(topo_dir, 'GH_tiles_2021_2s.asc')
        topofiles.append([3, topo_file])

        # 2 arc sec topo:
        #topo_file = os.path.join(topo_dir, 'astoriaS_2s_mhw.asc')
        topo_file = os.path.join(topo_dir, 'astoria_2s_mhw.asc')
        topofiles.append([3, topo_file])

    if 1:
        # 1/3 arc sec topo:
        topo_file = os.path.join(topo_dir,\
                                 'GH_tiles_2021_13s.asc')
        topofiles.append([3, topo_file])

    # ---------------
    # DTOPO:
    # ---------------
    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [dtopotype, fname]

    dtopo_data.dtopofiles = [[3, dtopofile]]
    dtopo_data.dt_max_dtopo = 10.


    # ---------------
    # qinit:
    # ---------------
    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    rundata.qinit_data.variable_eta_init = True  # for subsidence


    # ---------------
    # Force dry:
    # ---------------
    if 0:
        #None for this project
        namey = 'force_dry_init.data'
        force_dry_fname = os.path.join(input_dir, namey)
        force_dry = ForceDry()
        force_dry.tend = 1e9
        force_dry.fname = force_dry_fname
        rundata.qinit_data.force_dry_list.append(force_dry)


    # ---------------
    # REGIONS:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    flagregions = rundata.flagregiondata.flagregions  # initialized to []

    #Level 1 is 4', Level 2 is 2', Level 3 is 24", Level 4 is 12", Level 5 is 6",
    #Level 6 is 3", Level 7 is 1", Level 8 is 1/3", Level 9 is 1/9"
    # dx = dy = 4', 2', 24", 12", 6", 3", 1", 1/3", 1/9"
    #refinement_ratios = [2,5,2,2,2,3,3,3]


    # Computational domain Variable Region - 4 minutes to 24":
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
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

    # Region12sec - fixed 12 sec:
    # Level 4 is 12 sec
    # Note that this is a rectangle specified in the new way
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_12sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-126.6,-123.6,46.27,47.68]
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 24", 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51b_24sec_12sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    #flagregion.t2 = tmax_dtopo_region
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

    if 1:
        ### Will use this for the inundation runs. (6" slider window)
        # Rectangular region that encompasses gauges 94-137, offshore OSVES
        # or gauges 98-143, offshore Westport, or gauges 166-208 offshore Seaside.
        # Make this region 6" for all time to check if gauge values change
        # This rectangle goes out to the Ruled Rectangle (without the b) above.
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_6sec'
        flagregion.minlevel = 5
        flagregion.maxlevel = 5
        flagregion.t1 = 0.0
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle for now

        ## for Ocean Shores, encompasses gauges 94 to 137
        #latitudes below are for 93 in north to 138 in south
        gauges_region = [-125.9,-124.1,46.66,47.287]

        ## for Westport, encompasses gauges 98 to 143
        #latitudes below are for 97 in north to 144 in south
        #gauges_region = [-125.9,-124.1,46.587,47.227]

        ## for Hoquiam, encompasses gauges 94 to 143
        #latitudes below are for 93 in north to 144 in south
        #gauges_region = [-125.9,-123.82,46.6,47.227]

        ## for Seaside, encompasses gauges 166 to 208
        #gauges_region = [-125.65,-123.76,45.69,46.43]

        flagregion.spatial_region = gauges_region
        flagregions.append(flagregion)

    ## These are for WestportN:

    if 1:
        # 3sec region
        # Level 6 is 3 sec
        # Note that this is a rectangle specified in the new way
        # (other regions below will force/allow more refinement)
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_3sec_WestportN'
        flagregion.minlevel = 6
        flagregion.maxlevel = 6
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.245,-123.82, 46.83,47.06]
        flagregions.append(flagregion)

        # 1sec region
        # Level 7 is 1 sec
        # Note that this is a rectangle specified in the new way
        # (other regions below will force/allow more refinement)
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_1sec_WestportN'
        flagregion.minlevel = 7
        flagregion.maxlevel = 7
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.2,-124.03,46.85,46.95]
        flagregions.append(flagregion)

        # 1/3sec region
        # Level 8 is 1/3 sec
        # Note that this is a rectangle specified in the new way
        # (other regions below will force/allow more refinement)
        # Note that this region contains the fgmax_extent, but is larger
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_13sec_WestportN'
        flagregion.minlevel = 8
        flagregion.maxlevel = 8
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.18,-124.07,46.875,46.925]
        flagregions.append(flagregion)

    # ---------------
    # GAUGES:
    # ---------------

    if 0:
        rundata.gaugedata.gauges = []
        gauge_numbers = gauge_x.keys()
        for k in gauge_numbers:
            if gauge_flag[k]:
               rundata.gaugedata.gauges.append([k,gauge_x[k],gauge_y[k],0,1e9])

    # ---------------
    # FIXED GRIDS:
    # == setfixedgrids.data values ==
    # ---------------
    # Not in Geoclaw 10, taking out
    #fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # -----------------------------
    # FGMAX GRIDS:
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

    fgmax_extent=[-124.17, -124.09, 46.88, 46.92]

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
        fg.tstart_max =  t0 + 120   # when to start monitoring max values
        fg.tend_max = 1.e10         # when to stop monitoring max values
        fg.dt_check = 5.           # target time (sec) increment between updating
                                    # max values
                                    # which levels to monitor max on
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-2      # tolerance for flagging arrival

        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data

    # ======================================
    if 1:
        ###  HERE IS THE FGOUT STUFF TO EDIT (if plot_eta is set to True in
        ###  compare_gauge_max_withasce.py postprocessing, it uses fgout)
        # == fgout_grids.data values ==
        # NEW IN v5.9.0
        # Set rundata.fgout_data.fgout_grids to be a list of
        # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:

        #This is the region where 1/3 computation was
        #flagregion.spatial_region = [-124.185,-124.1325,46.9725,46.99]
        one_sixth = 1.0/(6.0*3600.)

        #fgout_extent = fgmax_extent

        fgout_extent = [-124.14, -124.07, 46.88, 46.92]

        fgout_grids = rundata.fgout_data.fgout_grids  #empty list initially
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 1
        fgout.point_style = 2  #a 2d grid of points
        fgout.output_format = 'binary32'
        dt = 10.
        dx_fine = 1./(3*3600.)  # grid resolution at finest level
        # NEED TO FIX:
        fgout.x1 = fgout_extent[0] #+ dx_fine/2.
        fgout.x2 = fgout_extent[1] #- dx_fine/2.
        fgout.y1 = fgout_extent[2] #+ dx_fine/2.
        fgout.y2 = fgout_extent[3] #- dx_fine/2.
        fgout.tstart = 1800.
        fgout.tend = 1.5*3600
        fgout.nout = int(np.round((fgout.tend - fgout.tstart)/dt)) + 1
        fgout.nx = int(np.round((fgout.x2-fgout.x1)/dx_fine))
        fgout.ny = int(np.round((fgout.y2-fgout.y1)/dx_fine))
        fgout_grids.append(fgout)    # written to fgout_grids.data



    # ---------------
    # DEVELOPERS:
    # ---------------
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
    amrdata.tprint = False       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    from clawpack.geoclaw import kmltools
    rundata = setrun(*sys.argv[1:])
    rundata.write()

    # To create kml files of inputs:
    #kmltools.make_input_data_kmls(rundata)
