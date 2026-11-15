"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""
import sys
import os
import numpy as np
import glob

# if some values in .data files show up as e.g. np.float64(3600.0)
# this will restore old behavior and just print 3600.0:
#np.set_printoptions(legacy="1.25")
# fixed in clawpack versions >= v5.13.0

from clawpack.geoclaw import topotools, dtopotools
from clawpack.geoclaw import fgmax_tools, fgout_tools

from clawpack.amrclaw.data import FlagRegion

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

try:
    root_dir = os.environ['GTT']
except:
    raise Exception("*** Must first set GTT enviornment variable")
    raise Exception("*** Perhaps to point to geoclaw_tsunami_tutorial/GTT")

topodir =  root_dir + '/topo/topofiles'
dtopodir = root_dir + '/dtopo/dtopofiles'
RRdir = root_dir + '/topo/regions'


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
    # Standard Clawpack parameters to be written to claw.data, amr.data:
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim      # should always be 2 in GeoClaw

    # Lower and upper edge of computational domain:

    domain_extent = (-137, -123, 38, 54)
    x_west, x_east, y_south, y_north = domain_extent  # unpack

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = x_west
    clawdata.upper[0] = x_east

    clawdata.lower[1] = y_south
    clawdata.upper[1] = y_north

    # Number of grid cells: Coarsest grid is 4 arcminutes (1/15 degree)
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

    clawdata.t0 = 0.  # Start time in seconds

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

    clawdata.output_style = 1

    if 0:
        # test timing of read
        clawdata.num_output_times = 0
        clawdata.tfinal = 1.0
        clawdata.output_t0 = False  # output at initial (or restart) time?


    elif clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 6
        clawdata.tfinal = 2.0*3600.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        # e.g. at 0 and 60 seconds and then less frequently:
        clawdata.output_times = [0.,60., 1800., 3600., 5400.]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True


    clawdata.output_format = 'binary32'      # 'ascii', 'binary', or 'binary32'
    clawdata.output_q_components = 'all'     # output all 3 components h,hu,hv


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
    clawdata.cfl_desired = 0.85
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

    clawdata.checkpt_style = 1

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = 1800.*np.arange(1,6,1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # maximum size of each grid patch (in each direction):
    amrdata.max1d = 60

    # initial size of work array for AMR patches:
    amrdata.memsize = 10000000

    # max number of refinement levels:
    amrdata.amr_levels_max = 6

    # List of refinement ratios at each level (length at least mxnest-1)

    # dx = dy = 4', 1', 30", 6", 2", 1", 1/3"":
    #refinement_ratios = [4,2,5,3,2,3]

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


    # ---------------
    # REGIONS:
    # ---------------

    flagregions = rundata.flagregiondata.flagregions  # empty list initially

    # grid resolutions:
    # dx = dy = 4', 2', 24", 12", 6", 3", 1", 1/3", 1/9"

    # Computational domain Variable Region:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1
    flagregion.maxlevel = 2
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    x1,y1 = clawdata.lower
    x2,y2 = clawdata.upper
    # Domain is [x1,x2,y1,y2], and add a buffer around it (not really needed):
    flagregion.spatial_region = [x1-0.2, x2+0.2, y1-0.2, y2+0.2]
    flagregions.append(flagregion)

    # Ruled rectangles covering dtopo region, variable 24" to 12"

    # Continential shelf extended to cover dtopo
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

    # Continential shelf extended to cover dtopo
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

    # Fixed 12 sec region S. WA coast
    extent_12sec = [-126.6,-123.6,46.27,47.68]
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_12sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = extent_12sec
    flagregions.append(flagregion)

    # Fixed 6" region offshore Grays Harbor
    extent_6sec = [-125.9,-124.1,46.587,47.227]
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_6sec'
    flagregion.minlevel = 5
    flagregion.maxlevel = 5
    flagregion.t1 = 0.0
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle for now
    flagregion.spatial_region = extent_6sec
    flagregions.append(flagregion)

    ## These are for WestportCoast:
    # [-124.165, -124.06, 46.75, 46.92]

    if 1:
        # 3sec region
        extent_3sec = [-124.245,-123.82, 46.7,47.06]
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_3sec_WestportN'
        flagregion.minlevel = 6
        flagregion.maxlevel = 6
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = extent_3sec
        flagregions.append(flagregion)

    if 1:
        # 1sec region
        extent_1sec = [-124.2,-123.97,46.7,46.95]
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_1sec_WestportCoast'
        flagregion.minlevel = 7
        flagregion.maxlevel = 7
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = extent_1sec
        flagregions.append(flagregion)

        # 1/3sec region
        extent_13sec = [-124.18,-124.06,46.75,46.925]
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_13sec_WestportCoast'
        flagregion.minlevel = 8
        flagregion.maxlevel = 8
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = extent_13sec
        flagregions.append(flagregion)

    # ------------------------------
    # Set data specific to GeoClaw:
    # ------------------------------

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.1


    geo_data = rundata.geo_data

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
    geo_data.speed_limit = 20.  # limit on speed sqrt(u**2 + v**2)


    # ---------------
    # TOPO:
    # ---------------
    # == topo.data values ==
    topo_data = rundata.topo_data

    topofiles = topo_data.topofiles   # empty list initially
    # for topography, append tuples/lists of the form:
    #    [topotype, fname]

    #topodir = '../../topo/topofiles'  # path to topofiles used below

    if 1:
        # 15-sec topo:
        topo = topotools.Topography()
        topo.path = '/Users/rjl/topo/etopo22_mosaic/etopo_2022_NW_Pacific.nc'
        topo.topo_type = 4  # netCDF
        topo.crop_extent = domain_extent  # computational domain
        topo.buffer = 1  # buffer to make sure domain is covered
        topo.coarsen = 2 # coarsen from 15 arcsec to 30 arcsec
        topo_data.topofiles.append(topo)

    # 3" Coastal Relief Model (CRM vol 8)
    topo = topotools.Topography()
    topo.path = '/Users/rjl/topo/crm/crm_vol8_2025.nc'
    topo.topo_type = 4  # netCDF
    topo.crop_extent = extent_6sec
    topo.buffer = 10  # buffer to make sure refined region is covered
    topo.coarsen = 2 # coarsen from 3 arcsec to 6 arcsec
    topo_data.topofiles.append(topo)

    topo = topotools.Topography()
    topo.path = '/Users/rjl/topo/crm/crm_vol8_2025.nc'
    topo.topo_type = 4  # netCDF
    topo.crop_extent = extent_3sec
    topo.buffer = 10  # buffer to make sure refined region is covered
    topo.coarsen = 1
    topo_data.topofiles.append(topo)


    # 1/9" CUDEM times: find the ones that intersect the fine-grid region:
            
    def intersection(topo, extent):
        topo.read_header()
        x1t,x2t,y1t,y2t = topo._extent
        x1e,x2e,y1e,y2e = extent

        x1 = max(x1t, x1e)
        x2 = min(x2t, x2e)
        dx_int = x2 - x1

        y1 = max(y1t, y1e)
        y2 = min(y2t, y2e)
        dy_int = y2 - y1

        if dx_int >= topo._delta[0] and dy_int >= topo._delta[1]:
            return [x1,x2,y1,y2]
        else:
            return None

        
    if 1:
        crop_extent = extent_1sec
        tiles = glob.glob('/Users/rjl/topo/sowa_mhw_2025/nc_tiles_topotools/*nc')
        print(f'+++ Found {len(tiles)} 1/9" tiles to check')

        # select only the tiles that overlap with crop_extent:
        for tile_path in tiles:
            topo = topotools.Topography()
            topo.path = tile_path
            topo.topo_type = 4  # netCDF
            intersect_extent = intersection(topo, crop_extent)
            if intersect_extent is not None:
                # intersection of crop_extent with topo is non-empty
                print(f'+++ adding {tile_path}:\n    {intersect_extent}')
                topo.crop_extent = intersect_extent
                topo.coarsen = 9 # coarsen from 1/9 arcsec to 1 arcsec
                topo.x_align = -124.
                topo.y_align = 47.
                topo_data.topofiles.append(topo)




    # ---------------
    # DTOPO:
    # ---------------
    # == setdtopo.data values ==
    # for moving topography, append lists of the form  [dtopo_type, fname]
    # to the initially empty list rundata.dtopo_data.dtopofiles:
    dtopo_data = rundata.dtopo_data

    if 0:
        #dtopodir = '../../dtopo/dtopofiles'  # path to dtopofile used below
        dtopofile = os.path.join(dtopodir, 'ASCE_SIFT_Region2.dtt3')
        dtopo_data.dtopofiles.append([3, dtopofile])
        dtopo_data.dt_max_dtopo = 0.2  # max timestep (sec) while topo is changing

    else:
        #dtopofile = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/testing_netcdf/BL16D.nc'
        dtopofile = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/dtopo30sec_nosubevents_kinokada/dtopofiles_nc/BL16D.nc'
        #dtopo_data.dtopofiles.append([4, dtopofile])

        #dtopofile = '/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_groundmotions/dtopo30sec/dtopofiles/BL16D.dtt3'
        #dtopo_data.dtopofiles.append([3, dtopofile])

        dtopo = dtopotools.DTopography()
        dtopo.path = dtopofile
        dtopo.dtopo_type = 4

        dtopo_data.dtopofiles.append(dtopo)

        dtopo_data.dt_max_dtopo = 10  # max timestep (sec) while topo is changing


    # ---------------
    # qinit:
    # ---------------
    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    rundata.qinit_data.variable_eta_init = True  # for subsidence/uplift


    # ---------------
    # Force dry:
    # ---------------
    if 0:
        # None for this project
        name = 'force_dry_init.data'
        force_dry_fname = os.path.join(input_dir, name)
        force_dry = ForceDry()
        force_dry.tend = 1e9
        force_dry.fname = force_dry_fname
        rundata.qinit_data.force_dry_list.append(force_dry)



    # ---------------
    # GAUGES:
    # ---------------

    gauges = rundata.gaugedata.gauges   # empty list initially
    # for gauges append tuples/lists of the form
    #   [gaugeno, x, y, t1, t2]

    gauges.append([101, -124.1895833, 47.1162500, 0, 1e9]) # slightly offshore
    gauges.append([102, -124.1804167, 47.1162500, 0, 1e9]) # onshore
    gauges.append([103, -124.1704167, 47.1162500, 0, 1e9]) # in river

    rundata.gaugedata.file_format = 'ascii'  # often use 'binary32'
    #rundata.gaugedata.min_time_increment = 5 # minimum seconds between outputs


    # -----------------------------
    # FGMAX GRIDS:
    # ------------------------------
    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2
    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list initially

    fgmax_extent = [-124.165, -124.06, 46.75, 46.92]

    if 1:
        # Points on a uniform 2d grid:
        dx_fine = 3./(3600.)  # grid resolution at finest level
        fg = fgmax_tools.FGmaxGrid()
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = fgmax_extent[0] #+ dx_fine/2.
        fg.x2 = fgmax_extent[1] #- dx_fine/2.
        fg.y1 = fgmax_extent[2] #+ dx_fine/2.
        fg.y2 = fgmax_extent[3] #- dx_fine/2.
        fg.dx = dx_fine
        fg.tstart_max =  120   # when to start monitoring max values
        fg.tend_max = 1.e10         # when to stop monitoring max values
        fg.dt_check = 5.           # target time (sec) increment between updating
                                    # max values
                                    # which levels to monitor max on
        fg.min_level_check = amrdata.amr_levels_max
        fg.arrival_tol = 1.e-2      # tolerance for flagging arrival

        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data


    # -----------------------------
    # FGOUT GRIDS:
    # ------------------------------

    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    if 1:
        #fgout_extent = fgmax_extent

        fgout_extent = [-124.165, -123.97, 46.73, 46.92]

        fgout_grids = rundata.fgout_data.fgout_grids  #empty list initially
        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 1
        fgout.point_style = 2  #a 2d grid of points
        fgout.output_format = 'binary32'
        dt = 30.
        dx_fine = 3./(3600.) 
        # NEED TO FIX:
        fgout.x1 = fgout_extent[0] #+ dx_fine/2.
        fgout.x2 = fgout_extent[1] #- dx_fine/2.
        fgout.y1 = fgout_extent[2] #+ dx_fine/2.
        fgout.y2 = fgout_extent[3] #- dx_fine/2.
        fgout.tstart = 1800.
        fgout.tend = clawdata.tfinal
        fgout.nout = int(np.round((fgout.tend - fgout.tstart)/dt)) + 1
        fgout.nx = int(np.round((fgout.x2-fgout.x1)/dx_fine))
        fgout.ny = int(np.round((fgout.y2-fgout.y1)/dx_fine))
        fgout_grids.append(fgout)    # written to fgout_grids.data

    # ---------------
    # FOR DEVELOPERS:
    # ---------------

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
    rundata = setrun(*sys.argv[1:])
    rundata.write()

    # To create kml files of inputs:
    from clawpack.geoclaw import kmltools
    kmltools.make_input_data_kmls(rundata, dtopo_contours=True)
