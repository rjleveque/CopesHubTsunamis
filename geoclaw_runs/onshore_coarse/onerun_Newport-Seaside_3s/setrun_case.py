"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

Note: this version has case as a parameter to setrun, a dictionary used to pass
in values that vary from case to case for doing a parameter sweep.
"""

import os
import numpy as np
from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgout_tools, fgmax_tools

use_bouss_version = False  # True ==> requires code compiled with PETSc etc.
                           # might still solve SWE if bouss_equations = 0 below

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

root_dir = os.environ['CHT']
input_files_dir = root_dir + '/input_files/'
print('setting input_files_dir = ',input_files_dir)

rundir = os.getcwd()
print('rundir = %s' % rundir)

if 1:
    # for new ground motions:
    event = os.path.split(os.getcwd())[-1]
    instant = ('instant' in rundir)  # is this instantaneous uplift?
else:
    # for T-shirt event, hardwire here:
    event = 'SM1'
    instant = True

print('The event is  %s, instant = %s' % (event,instant))

#Note, we are forcing the outer RR to stay on, not just because of the
#dtopo, but because we want the wave to propogate to the inner RR region
#before letting go of the outer RR.

if instant:
    tmax_dtopo_region = 15*60.  # force fine grids up to this time
                                # 15 minutes for wave to propagate
                                # into inner RR before letting go of
                                # outer RR. Even in South, waves start
                                # to propogate immediately
    print ('assuming this is an instantaneous rupture')
    print (' ')

else:
    tmax_dtopo_region = 20*60.  # force fine grids up to this time
                                # 20 minutes for wave to propagate
                                # into inner RR in the South
    print ('assuming this is a kinematic rupture')
    print (' ')


topodir = root_dir + '/topo/topofiles'                        #  on Loyces laptop
#dtopodir = root_dir + '/dtopo/CSZ_groundmotions/dtopofiles'   # Not used it seems

# for hyak cluster:
topodir = topodir.replace('/mmfs1/home', '/gscratch/tsunami')
#dtopodir = dtopodir.replace('/mmfs1/home', '/gscratch/tsunami')


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

    # The values below are expected to be in case dictionary,
    # and may vary from case to case:

    try:
        dtopofiles = case['dtopofiles']
    except:
        # for make_kmls
        dtopofiles = [[3,'/Users/rjl/git/CopesHubTsunamis/dtopo/CSZ_Bandon_noext/CSZ_M1_noext.tt3']]

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
    one_sixth = 1.0/(3600.*6)
    one_third = 1.0/(3600.*3)


    #Shift the domain so the center of cells align with the etopo22 data
    clawdata.lower[0] = -135.       # west longitude
    clawdata.upper[0] = -122.       # east longitude
    clawdata.lower[1] = 38.5        # south latitude
    clawdata.upper[1] = 53.5        # north latitude

    clawdata.num_cells[0] = 13*15
    clawdata.num_cells[1] = 15*15


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
        ## ADJUST:
        clawdata.num_output_times = 6    # no frame output, only fgmax
        clawdata.tfinal = 3600.
        clawdata.output_t0 = False        # output at initial (or restart) time?

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

    if instant: #instantaneous rupture
        clawdata.dt_initial = 0.2
    else:       # kinematic rupture
        clawdata.dt_initial = 10.0

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
    #amrdata.amr_levels_max = 4
    amrdata.amr_levels_max = 6

    # List of refinement ratios at each level (length at least mxnest-1)

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
    #wont need this one since 15-second below covers whole domain
    #topofiles.append([3, topodir + '/etopo22_1min_-163_-122_38_63.asc'])

    # 15-second etopo22:
    #This will be sufficient for runs up to 5" for initial runs before
    #inundation runs.
    topofiles.append([3, topodir + '/etopo22_15s_-137_-121_37_55.asc'])


    if 1:
        # 3-second topo:
        topofiles.append([3, topodir + '/crm_vol8_3sec_cropped_44-47.asc'])


    # 2sec topo files coarsened from 1/9" navd88 tiles and adjusted
    # to approximate mhw by subtracting 2.1m from Z:

    topo2sec_files = \
        ['ncei19_n45x00_w124x00_amhw_2021_2sec.asc',
         'ncei19_n44x75_w124x00_amhw_2021_2sec.asc',
         'ncei19_n45x25_w124x00_px_amhw_2019_2sec.asc',
         'ncei19_n46x00_w124x00_px_amhw_2019_2sec.asc',
         'ncei19_n46x25_w124x00_px_amhw_2019_2sec.asc',
         'ncei19_n45x75_w124x25_px_amhw_2019_2sec.asc',
         'ncei19_n44x75_w124x25_amhw_2021_2sec.asc',
         'ncei19_n45x50_w124x25_px_amhw_2019_2sec.asc',
         'ncei19_n45x00_w124x25_amhw_2021_2sec.asc',
         'ncei19_n45x25_w124x25_px_amhw_2019_2sec.asc',
         'ncei19_n46x25_w124x25_px_amhw_2019_2sec.asc',
         'ncei19_n46x00_w124x25_px_amhw_2019_2sec.asc',
         'ncei19_n45x50_w124x00_px_amhw_2019_2sec.asc',
         'ncei19_n45x75_w124x00_px_amhw_2019_2sec.asc']
    for fname in topo2sec_files:
        topofiles.append([3, os.path.join(topodir, fname)])


    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]

    dtopo_dir = os.path.join(root_dir, 'dtopo/CSZ_groundmotions/dtopofiles')
    dtopo_data.dtopofiles = dtopofiles # from case dictionary


    if instant: #instantaneous rupture
        dtopo_data.dt_max_dtopo = .2
    else:       #dynamic rupture
        dtopo_data.dt_max_dtopo = 10.0


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

    # Computational domain Variable Region - 4min, 2min to 24 sec:
    # Level 3 below is 24 sec
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

    # Continential shelf extended to cover dtopo, 24"  #12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51b_24sec'
    flagregion.minlevel = 3  #4
    flagregion.maxlevel = 3  #4
    flagregion.t1 = 0.
    flagregion.t2 = tmax_dtopo_region
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51b.data')
    flagregions.append(flagregion)

    # Continential shelf extended to cover dtopo, 24"  #12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46b_24sec'
    flagregion.minlevel = 3  #4
    flagregion.maxlevel = 3  #4
    flagregion.t1 = 0.
    flagregion.t2 = tmax_dtopo_region
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46b.data')
    flagregions.append(flagregion)

    # New Continential shelf around onshore region, 12"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_near_onshore_12sec'
    flagregion.minlevel = 3  #4
    flagregion.maxlevel = 3  #4
    flagregion.t1 = 0.
    flagregion.t2 = tmax_dtopo_region
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-125.9, -123.8, 44.0, 48]
    flagregions.append(flagregion)

    # coastal region for fgmax:
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_3sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 6
    flagregion.t1 = 0.0
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle for now

    ## Newport to Seaside
    onshore_extent = [-124.1,-123.8,44.5,46.05] # also used for fgmax grid
    coast_near_onshore_region = [-124.15,-123.75,44.45,46.1]

    flagregion.spatial_region = coast_near_onshore_region
    flagregions.append(flagregion)


    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    #rundata.gaugedata.gauges.append([1,-122.4, 47.781, 0., 1.e10])

    rundata.gaugedata.file_format = 'binary32'

    if 0:
        asce_gauges_file = root_dir + '/gauges/offshore_gauges/asce_values.txt'
        asce_gauges = np.loadtxt(asce_gauges_file, skiprows=1)

        #Use all the ASCE gauges for this job run
        for k in range(0,len(asce_gauges)):
            gaugeno = int(asce_gauges[k,0])
            gx = float(asce_gauges[k,1])
            gy = float(asce_gauges[k,2])
            if gy >= 40:
                rundata.gaugedata.gauges.append([gaugeno, gx, gy, 0, 1e9])
        VI_gauges_file = root_dir + '/gauges/offshore_gauges/gaugesVI.txt'
        VI_gauges = np.loadtxt(VI_gauges_file, skiprows=1)
        for k in range(0,len(VI_gauges),1):
            gaugeno = int(VI_gauges[k,0])
            gx = float(VI_gauges[k,1])
            gy = float(VI_gauges[k,2])
            rundata.gaugedata.gauges.append([gaugeno, gx, gy, 0, 1e9])


    if 1:
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
        dx_fine = 3./3600.  # grid resolution at finest level for the big fgmax region


        #keep the domain shifted to align with the etopo22 data
        #This means that whole numbers like -130 and -123 below
        #will be at the center of 1/3" cells. They will also be
        #at centers of 30" cells (domain originally -135 to -122.
        #Likewise, 39.0 and 52.0 will be at centers of both 1/3"
        #and 30" cells (domain originally was 38.5 to 53.5)

        fg = fgmax_tools.FGmaxGrid()
        x1,x2,y1,y2 = onshore_extent
        fg.point_style = 2  # uniform rectangular x-y grid
        fg.x1 = x1 + dx_fine/2.0       # specified pts to align with FV cell centers
        fg.x2 = x2 - dx_fine/2.0       # since whole numbers are domain edges here
        fg.y1 = y1  + dx_fine/2.0
        fg.y2 = y2  - dx_fine/2.0
        fg.dx = dx_fine
        fg.dy = dx_fine
        fg.tstart_max =  0.       # when to start monitoring max values
        fg.tend_max = 1.e10       # when to stop monitoring max values
        fg.dt_check = 10.         # target time (sec) increment between updating
        #fg.min_level_check = 5   # monitor on finest level for a 5 level run
        fg.min_level_check = 3    # Coastal regions are 3,4 after tmax sec
        fg.arrival_tol = 1.e-1    # tolerance for flagging arrival
        fg.interp_method = 0      # 0 ==> pw const in cells, recommended
        fgmax_grids.append(fg)    # written to fgmax_grids.data

        #Add the fg max grid here at the finest level

    if 0:
        ###  HERE IS THE FGOUT STUFF TO EDIT (if plot_eta is set to True in
        ###  compare_gauge_max_withasce.py postprocessing, it uses fgout)
        # == fgout_grids.data values ==
        # NEW IN v5.9.0
        # Set rundata.fgout_data.fgout_grids to be a list of
        # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
        fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially
        dx_fgout = 1/60.  # degrees
        dy_fgout = dx_fgout
        dt_fgout = 30  # seconds

        fgout = fgout_tools.FGoutGrid()
        fgout.fgno = 1
        fgout.point_style = 2       # will specify a 2d grid of points
        fgout.output_format = 'binary32'
        #fgout.nx = 8*60
        #fgout.ny = 12*60

        fgout.x1 = -130.   # specify edges (fgout pts will be cell centers)
        fgout.x2 = -122.   # edge of a cell, edges are whole numbers in this domain
        fgout.y1 = 38.5    # edge of a cell
        fgout.y2 = 50.5    # edge of a cell
        fgout.nx = int(round((fgout.x2 - fgout.x1)/dx_fgout))
        fgout.ny = int(round((fgout.y2 - fgout.y1)/dy_fgout))
        fgout.tstart = 0.
        fgout.tend = 1*3600
        fgout.nout = int(round(((fgout.tend - fgout.tstart)/dt_fgout))) + 1
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
    rundata = setrun(*sys.argv[1:])
    rundata.write()

    if 1:
        from clawpack.geoclaw import kmltools
        sys.path.insert(0,root_dir+'/common_code')  # works for fgmax pt style 4
        import kmltools
        kmltools.make_input_data_kmls(rundata)
