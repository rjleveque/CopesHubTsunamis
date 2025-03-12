"""
Given a point (x_desired, y_desired) where you want a gauge, roughly.

Locate a point (x_center, y_center) at the center of a grid cell
at the resolution specified by (dx, dy), that is as close as possible to the
desired point.
"""

import numpy

def adjust(z_desired, z_edge, dz, verbose=False):
    """
    Given a desired location z (either x or y) for a gauge or other point
    of interest, adjust the location is it is offset (integer + 1/2)*dz
    from z_edge, an arbitrary edge location (e.g. an edge of the
    computational domain or any point offset integer*dz from the domain edge).
    This will put the new location in the center of a finite volume cell
    provided this point is in a patch with resolution dz.
    """
    # z represents either x or y
    i = numpy.round((z_desired-z_edge - 0.5*dz)/dz)
    z_centered = z_edge + (i+0.5)*dz
    if verbose:
        zshift = z_centered - z_desired
        zfrac = zshift / dz
        print("adjusted from %15.9f" % z_desired)
        print("           to %15.9f" % z_centered)
        print("   shifted by %15.9f = %.3f*dz" % (zshift,zfrac))
    return z_centered

def adjust_xy(x_desired, y_desired, x_edge, y_edge, dx, dy, verbose=False):
    assert len(x_desired) == len(y_desired), \
            '*** lengths of x_desired, y_desired do not match'
    x_centered = []
    y_centered = []
    for i in range(len(x_desired)):
        x = adjust(x_desired[i], x_edge, dx, verbose=verbose)
        y = adjust(y_desired[i], y_edge, dy, verbose=verbose)
        x_centered.append(x)
        y_centered.append(y)
        
    return numpy.array(x_centered), numpy.array(y_centered)
        
def test():
    dx = 1/(3*3600.)
    xd = [-122.01, -122.02]
    yd = [47.001, 47.002]
    x_edge = -123.
    y_edge = 47.
    print('Desired x = ',xd)
    print('Desired y = ',yd)
    xc,yc = adjust_xy(xd,yd,x_edge,y_edge,dx,dx,True)
    print('Centered x = ',xc)
    print('Offsets in x in units of 1/3 arcsec: ', (xc-x_edge)*3*3600)
    print('Centered y = ',yc)
    print('Offsets in y in units of 1/3 arcsec: ', (yc-y_edge)*3*3600)

if __name__ == '__main__':
    test()
