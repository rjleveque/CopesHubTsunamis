"""
Code to track debris paths as points.
Set up for SouthBeach area, with initial debris given by:
    tracer particles in the water,
    particles on land that have a grounding depth

"""
from pylab import *

def compute_debris_paths(db_times, qoi_fcns):

    """
    Initial debris locations are then set, along with parameters such
    as mass and grounding depth.

    Code from debris_tools is used to move the debris based on the velocities
    loaded from the fgout frames.

    Returns debris_paths, dbnos_water, dbnos_land
    """

    #import debris_tools
    from geoclaw_debris import debris_tools

    #from clawpack.geoclaw import fgout_tools
    from geoclaw_debris import fgout_tools

    from scipy.interpolate import interp1d



    # First time, to initialize particles
    t0 = db_times[0]

    # initial topography (used to determine land / water points):
    # assumed to be same as B in fgout frames (no subsidence)
    B0_fcn = lambda x,y: qoi_fcns['B'](x,y,t0)  # function of (x,y) alone


    # define beach region for placing initial particles only on beach:


    # Initialize debris_paths dictionary and set
    # initial debris particle locations (x0,y0) at time t0.
    # Require a list of dbnos and each array
    #     debris_paths[dbno]
    # in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

    debris_paths = {}
    dbnos_water = []
    dbnos_land = []
    dbnos = []

    # same for all debris particles:
    grounding_depth = {}
    grounding_speed = 0.
    drag_factor = {}
    #drag_factor = 10.

    # set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
    u0 = 0.
    v0 = 0.

    # rectangular region to place initial particles:
    x1,x2,y1,y2 = [-124.14, -124.07, 46.74, 46.92]
    Bshore = 0.  # elevation of shore defining land vs water

    dxd = 0.002  # spacing between particles in latitude
    dyd = dxd * cos(y1*pi/180)  # spacing in longitude
    xd = arange(x1, x2, dxd)
    yd = arange(y1, y2, dyd)
    mxd = len(xd)
    myd = len(yd)

    for i in range(len(xd)):
        x0 = xd[i]
        for j in range(len(yd)):
            y0 = yd[j]
            dbno = myd*i + j
            db = array([[t0, x0, y0, u0, v0]])
            debris_paths[dbno] = db

            if B0_fcn(x0,y0) < Bshore:
                dbnos_water.append(dbno)
                grounding_depth[dbno] = 0.
                drag_factor[dbno] = None
            else:
                dbnos_land.append(dbno)
                grounding_depth[dbno] = 0.5
                drag_factor[dbno] = None #0.1

    dbnos = dbnos_water + dbnos_land
    print(f'Created {len(dbnos)} initial debris particles')
    print(f'      {len(dbnos_water)} in the water and {len(dbnos_land)} onshore')


    # Compute debris path for each particle by using all the fgout frames
    # in the list fgframenos (first frame should be one used to set t0 above):

    debris_paths = debris_tools.make_debris_paths_fcnxyt(db_times,
                                            qoi_fcns, 
                                            debris_paths, dbnos, drag_factor, grounding_depth, grounding_speed)

    # done computing debris paths
    # ------------------------------

    return debris_paths, dbnos_water, dbnos_land
