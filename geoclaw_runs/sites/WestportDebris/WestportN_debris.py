from pylab import *

def compute_debris_paths(tfinal):

    """
    Load in fgout frames up to time tfinal from the tsunami simulation.

    Initial debris locations are given by building locations.

    Code from debris_tools is used to move the debris based on the velocities
    loaded from the fgout frames.
    """
    import debris_tools
    #from clawpack.geoclaw import debris_tools
    from clawpack.geoclaw import fgout_tools
    from scipy.interpolate import interp1d

    outdir = '_output'
    output_format = 'binary32'
    qmap = 'geoclaw'  # defines mapping for fgout variables

    fgno = 1
    fgout_grid2 = fgout_tools.FGoutGrid(fgno, outdir, output_format, qmap)
    fgout_grid2.read_fgout_grids_data(fgno)
    #fgout_grid2.read_fgout_grids_data_pre511(fgno)
    print('Looking for output in ',outdir)

    # use all frames up to time tfinal:
    fgframes2 = [n+1 for n in range(len(fgout_grid2.times)) \
                   if fgout_grid2.times[n] <= tfinal]

    #fgframes2 = range(1,182,3) # incomplete run


    fgframe2 = fgframes2[0] # initial frame for fgout grid 2
    fgout2 = fgout_grid2.read_frame(fgframe2)

    # Deterime time t of first fgout frame, to initialize particles
    t0 = fgout2.t
    fgout_extent = fgout2.extent_edges

    # initial topography (used to determine land / water points):
    # assumed to be same as B in fgout frames (no subsidence)
    B0_fcn = fgout_tools.make_fgout_fcn_xy(fgout2, 'B')

    # define beach region for placing initial particles only on beach:


    # Initialize debris_paths dictionary and set
    # initial debris particle locations (x0,y0) at time t0.
    # Require a list of dbnos and each array
    #     debris_paths[dbno]
    # in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

    debris_paths = {}
    #dbnos_water = []
    #dbnos_land = []
    dbnos = []

    # same for all debris particles:
    grounding_depth = 0.
    grounding_speed = 0.
    drag_factor = None
    #drag_factor = 10.

    # set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
    u0 = 0.
    v0 = 0.

    fname_buildings = '/Users/rjl/git/CopesHubTsunamis/gauges/GraysHarborPacificCountyBuildings/WestportN_buildings.csv'
    fid_building,x_building,y_building = loadtxt(fname_buildings, delimiter=',',
                                                 unpack=True)
    for k in range(len(fid_building)):
        dbno = fid_building[k]
        db = array([[t0, x_building[k], y_building[k], u0, v0]])
        debris_paths[dbno] = db

        #    if B0xy < 0:
        #        dbnos_water.append(dbno)
        #    elif B0xy > 0:
        #        dbnos_land.append(dbno)
        dbnos.append(dbno)

    #dbnos = dbnos_water + dbnos_land
    print('Created %i initial debris particles' % len(dbnos))
    #print('    %i offshore,  %i onshore' % (len(dbnos_water),len(dbnos_land)))

    # Compute debris path for each particle by using all the fgout frames
    # in the list fgframes (first frame should be one used to set t0 above):

    debris_paths = debris_tools.make_debris_paths(fgframes2, fgout_grid2,
                                debris_paths, dbnos, drag_factor, grounding_depth, grounding_speed)

    # done computing debris paths
    # ------------------------------
    return debris_paths, dbnos
