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

instant = ('instant' in rundir)  # is this instantaneous uplift?

if instant:
    tmax_dtopo_region = 10.  # force fine grids up to this time
else:
    #tmax_dtopo_region = 300.   # force fine grids up to this time
    tmax_dtopo_region = 15*60.  # force fine grids up to this time
                                # 15 minutes for wave to propagate


if '/projects' in rundir:
    topodir = '/projects/rale6846/topo/topofiles'  # on CU
    dtopodir = '/projects/rale6846/dtopo/dtopofiles'  # on CU
else:
    #topodir = '/Users/rjl/topo/topofiles'        # on Randys laptop
    #dtopodir = '/Users/rjl/B/dtopo/dtopofiles'   # on Randys laptop
    topodir = root_dir + '/topo/topofiles'             # on Loyces laptop
    dtopodir = root_dir + '/dtopo/CSZ_groundmotions'   # on Loyces laptop


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
    one_sixth = 1.0/(3600.*6)
    one_third = 1.0/(3600.*3)


    #Shift the domain so the center of cells align with the etopo22 data
    clawdata.lower[0] = -135. - one_sixth      # west longitude
    clawdata.upper[0] = -122. - one_sixth      # east longitude
    clawdata.lower[1] = 38.5  - one_sixth      # south latitude
    clawdata.upper[1] = 53.5  - one_sixth      # north latitude

    clawdata.num_cells[0] = 13*15
    clawdata.num_cells[1] = 15*15


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
        clawdata.num_output_times = 9    #output every 5 minutes
        #clawdata.tfinal = 2*3600.
        clawdata.tfinal = 2700.          #run for 45 minutes
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
    amrdata.amr_levels_max = 4
    #amrdata.amr_levels_max = 5

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
    refinement_ratios = [4,2,2,3,5,3,3]


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

    dtopo_data.dtopofiles = [[3, dtopodir + '/buried-locking-str10-deep.dtt3']]

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
    
    # Computational domain Variable Region - 4min, 1min to 30 sec:
    # Level 3 below is 30 sec 
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

    # Continential shelf extended to cover dtopo, 15"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51b_15sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = tmax_dtopo_region
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51b.data')
    flagregions.append(flagregion)
    
    # Continential shelf extended to cover dtopo, 15"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46b_15sec'
    flagregion.minlevel = 4
    flagregion.maxlevel = 4
    flagregion.t1 = 0.
    flagregion.t2 = tmax_dtopo_region
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46b.data')
    flagregions.append(flagregion)

    # Continential shelf Variable Region, 30" to 15"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_46_51_30to15sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = tmax_dtopo_region
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_46_51.data')
    flagregions.append(flagregion)
    
    # Continential shelf Variable Region, 30" to 15"
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Coast_40_46_30to15sec'
    flagregion.minlevel = 3
    flagregion.maxlevel = 4
    flagregion.t1 = tmax_dtopo_region
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(RRdir + \
            '/RuledRectangle_Coast_40_46.data')
    flagregions.append(flagregion)

    if 1: #For 15" run around region of interest
        # Rectangular region that encompasses gauges 94-137, offshore OSVES
        # or gauges 98-143, offshore Westport.
        # Make this region 15" for all time thinking gauge plots will be better
        # and tsunami propagating in our sliver of interest.
        # This rectangle goes out to the Ruled Rectangle (without the b) above.
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_15sec'
        flagregion.minlevel = 4
        flagregion.maxlevel = 4
        flagregion.t1 = tmax_dtopo_region
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle for now

        ## for Ocean Shores, encompasses gauges 94 to 137
        #latitudes below are for 93 in north to 138 in south
        #gauges_region = [-125.9,-124.1,46.66,47.287]

        ## for Westport, encompasses gauges 98 to 143
        #latitudes below are for 97 in north to 144 in south
        gauges_region = [-125.9,-124.1,46.587,47.227]

        flagregion.spatial_region = gauges_region
        flagregions.append(flagregion)


    if 0: #For the 5" runs
        # Rectangular region out to the b coastal region from the destination
        # area for the time the source is moving to get the peak at 5".
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_5sec_initially'
        flagregion.minlevel = 5
        flagregion.maxlevel = 5
        flagregion.t1 = 0.
        flagregion.t2 = tmax_dtopo_region
        flagregion.spatial_region_type = 1  # Rectangle for now

        ## for Ocean Shores, encompasses gauges 94 to 137
        #latitudes below are for 93 in north to 138 in south
        #source_region = [-126.58,-124.1,46.66,47.287]

        ## for Westport, encompasses gauges 98 to 143
        #latitudes below are for 97 in north to 144 in south
        source_region = [-126.58,-124.1,46.587,47.227]
        flagregion.spatial_region = source_region
        flagregions.append(flagregion)

        # Rectangular region that encompasses gauges 94-137, offshore OSVES
        # or gauges 98-143, offshore Westport.
        # Make this region 5" for all time to check if gauge values change
        # from when these gauges were in the 3,4 region above and had at most 15"
        # This rectangle goes out to the Ruled Rectangle (without the b) above.
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_5sec'
        flagregion.minlevel = 5
        flagregion.maxlevel = 5
        flagregion.t1 = tmax_dtopo_region
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle for now

        ## for Ocean Shores, encompasses gauges 94 to 137
        #latitudes below are for 93 in north to 138 in south
        #gauges_region = [-125.9,-124.1,46.66,47.287]

        ## for Westport, encompasses gauges 98 to 143
        #latitudes below are for 97 in north to 144 in south
        gauges_region = [-125.9,-124.1,46.587,47.227]

        flagregion.spatial_region = gauges_region
        flagregions.append(flagregion)

    if 0: #will need to fix the regions below for each collaboratory
          # dx = dy = 4', 1', 30", 15", 5", 1", 1/3", 1/9"
    
        # Westport Variable Region - 30sec to 15sec to 5sec:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Westport_30-15-5sec'
        flagregion.minlevel = 3
        flagregion.maxlevel = 5
        flagregion.t1 = 0.
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.29, -123.655, 46.325, 47.159]
        flagregions.append(flagregion)

        # Grays Harbor Region - allow  1" grids, require 15sec at least:
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'Region_Grays_15-5-1sec'
        flagregion.minlevel = 4
        flagregion.maxlevel = 6
        flagregion.t1 = tstart_finestgrid
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.199, -123.809, 46.8, 47.145]
        flagregions.append(flagregion)

        # fixedgrid Region - require  5" grids, allow 1 and  1/3":
        flagregion = FlagRegion(num_dim=2)
        flagregion.name = 'fixedgrid_5-1-13sec'
        flagregion.minlevel = 5
        flagregion.maxlevel = 7
        flagregion.t1 = tstart_finestgrid
        flagregion.t2 = 1e9
        flagregion.spatial_region_type = 1  # Rectangle
        flagregion.spatial_region = [-124.185,-123.935,46.785,46.935]
        flagregions.append(flagregion)

        if 0:
            #Will need a 7x7 perhaps
            # fixedgrid Region - require  1/3":
            flagregion = FlagRegion(num_dim=2)
            flagregion.name = 'fixedgrid_13sec'
            flagregion.minlevel = 7
            flagregion.maxlevel = 7
            flagregion.t1 = tstart_finestgrid
            flagregion.t2 = 1e9
            flagregion.spatial_region_type = 1  # Rectangle
            flagregion.spatial_region = [-124.185,-123.935,46.785,46.935]
            flagregions.append(flagregion)

        if 0:
            #Will need an 8x8 perhaps
            # fixedgrid Region - require  1/9":
            flagregion = FlagRegion(num_dim=2)
            flagregion.name = 'fixedgrid_19sec'
            flagregion.minlevel = 8
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

    asce_gagues_file = root_dir + '/info/asce_values.txt'
    asce_gauges = np.loadtxt(asce_gagues_file, skiprows=1)
    for k in range(0,len(asce_gauges),10):
        gaugeno = int(asce_gauges[k,0])
        gx = float(asce_gauges[k,1])
        gy = float(asce_gauges[k,2])
        if gy >= 40:
            rundata.gaugedata.gauges.append([gaugeno, gx, gy, 0, 1e9])
    VI_gagues_file = root_dir + '/info/gaugesVI.txt'
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
    dx_fine = 30./3600.  # grid resolution at finest level for the big fgmax region
                         # use finest of 30"


    #keep the domain shifted to align with the etopo22 data
    #This means that whole numbers like -130 and -123 below
    #will be at the center of 1/3" cells. They will also be
    #at centers of 30" cells (domain originally -135 to -122.
    #Likewise, 39.0 and 52.0 will be at centers of both 1/3"
    #and 30" cells (domain originally was 38.5 to 53.5)

    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid
    fg.x1 = -130.      # specify pts to align with FV cell centers
    fg.x2 = -123. 
    fg.y1 = 39.0 
    fg.y2 = 52.0 
    fg.dx = dx_fine
    fg.dy = dx_fine
    fg.tstart_max =  0.       # when to start monitoring max values 
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 10.         # target time (sec) increment between updating
    #fg.min_level_check = 5    # monitor on finest level
    fg.min_level_check = 3     # Coastal regions are 3,4 after 300 sec, and 4,4 for 300 sec 
    fg.arrival_tol = 1.e-1    # tolerance for flagging arrival
    fg.interp_method = 0      # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)    # written to fgmax_grids.data

    #Add the fg max grid here at the finest level

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
    fgout.x1 = -130. - one_sixth  # specify edges (fgout pts will be cell centers)
    fgout.x2 = -122. - one_sixth  # edge of a cell
    fgout.y1 = 38.5  - one_sixth  # edge of a cell
    fgout.y2 = 50.5  - one_sixth  # edge of a cell
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
