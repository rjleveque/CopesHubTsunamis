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
np.set_printoptions(legacy="1.25")
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
        print('**** no dtopofile specified ****')

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


    one_sixth = 1.0/(6.0*3600.)
    #one_sixth = 0.  # FIX
    clawdata.lower[0] = -135. - one_sixth      # west longitude
    clawdata.upper[0] = -123.6- one_sixth      # east longitude
    clawdata.lower[1] = 38.5  - one_sixth      # south latitude
    clawdata.upper[1] = 51.5  - one_sixth      # north latitude

    clawdata.num_cells[0] = 171 # 11.4*15
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
        # Here we run for 120 seconds and produce NO time frame output:
        clawdata.num_output_times = 0
        clawdata.tfinal = 120.
        clawdata.output_t0 = False

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

    clawdata.checkpt_style = 0

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
    #amrdata.amr_levels_max = 6   # for 1/3" resolution
    #amrdata.amr_levels_max = 4    # for 3"
    amrdata.amr_levels_max = 1  # for short test

    # List of refinement ratios at each level (length at least mxnest-1)


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
    refinement_data.wave_tolerance = 0.05

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]

    # 15-second etopo22:
    topofiles.append([3, f'{topo_dir}/etopo22_15s_-137_-121_37_55.asc'])


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

    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    # one sample gauge near hypocenter of BL10D and BL10M:
    gauges.append([101, -126.3, 48., 0., 1e9])


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
    from clawpack.geoclaw import kmltools

    # this main program is only used if you run this as a script via:
    #   python setrun_case.py
    # or via
    #   make data
    # not when runclaw_makeplots_dtopos.py is used for multiple cases

    # to make data without any dtopofile, to create kmls for other inputs,
    # set case={}:
    rundata = setrun('geoclaw', case={})
    rundata.write()

    # To create kml files of inputs:
    kmltools.make_input_data_kmls(rundata)
