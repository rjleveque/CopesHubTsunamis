"""
Take a set of gauges specified as in setrun.py and adjust them to be
centered in grid cells that are dx by dy and that are an integer number
of dx or dy away from x_edge and y_edge, respectively.

Then print the new gauge locations out in a way that can be pasted back
into setrun.py.

See  https://www.clawpack.org/nearshore_interp.html
for details on why this might be needed.

See also the documentation in the module center_points.py used below,
which can be used more generally for centering arbitrary points.

Typically (dx,dy) should be the finest refinement level expected at the
gauge location over the time period when the gauge output will be used.
Note that (dx,dy) may need to be chosen differently for different gauges.

Note that the domain need not be refined to this level everywhere in order
for grids at resolution (dx,dy) to be offset from the edges of the domain
by integer multiples of (dx,dy).

The initial "desired" points may be shifted by as much as (dx/2, dy/2) to
achieve this centering.
"""

import numpy as np
import center_points

# original "desired" locations:
gauges = []
if 1:
    #Quileute Marina
    gauges.append([1, -124.636793, 47.911337, 0, 1e+09])

    #Quileute Human Services
    gauges.append([2, -124.635659, 47.911024, 0, 1e+09])

    #Quileute Tribal Center
    gauges.append([3, -124.636372, 47.909623, 0, 1e+09])

    #USGS Station 
    gauges.append([4, -124.634667, 47.913157, 0, 1e+09])

    #Rialto Beach 
    gauges.append([5, -124.638281, 47.917166, 0, 1e+09])

    #Salty Heifers 
    gauges.append([6, -124.638966, 47.908110, 0, 1e+09])

    #Quileute Oceanside Resort
    gauges.append([7, -124.632049, 47.906042, 0, 1e+09])

    #Quileute Tribal School
    gauges.append([8, -124.610552, 47.891274, 0, 1e+09])

if 0:
    #Bureau Indian Affairs
    gauges.append([1, -124.419231, 47.743728, 0, 1e+09])

    #Hoh River mouth
    gauges.append([2, -124.433699, 47.748427, 0, 1e+09])

    #Lower Hoh Rd - River entrance
    gauges.append([3, -124.432322, 47.748407, 0, 1e+09])

    #2455 Lower Hoh Rd
    gauges.append([4, -124.422210, 47.745031, 0, 1e+09])

    #Lower Hoh Rd-Chalatt Loop
    gauges.append([5, -124.416469, 47.742495, 0, 1e+09])

    #Chalaat Loop
    gauges.append([6, -124.413932, 47.742409, 0, 1e+09])

    #Chief Klia Wellness Center
    gauges.append([7, -124.413180, 47.740973, 0, 1e+09])

    #Oil City Trailhead
    gauges.append([8, -124.408580, 47.749361, 0, 1e+09])

    #9772 Oil City Rd
    gauges.append([9, -124.404040, 47.749936, 0, 1e+09])

    #Oil City Road - Fossil Creek
    gauges.append([10, -124.405132, 47.749701, 0, 1e+09])

    #National Playground
    gauges.append([11, -124.415120, 47.741840, 0, 1e+09])

    #Pacific NW Trail
    gauges.append([12, -124.424083, 47.749583, 0, 1e+09])

    #Lower Hoh Rd East
    gauges.append([13, -124.403684, 47.740672, 0, 1e+09])

    #West Shore
    gauges.append([14, -124.427831, 47.742840, 0, 1e+09])

    #US 101 at Ruby Beach
    gauges.append([15, -124.4111, 47.71, 0, 1e+09])

    #US 101 at Kalaloch
    gauges.append([16, -124.37158, 47.6071, 0, 1e+09])


if 0:  #Desired at the moment, might need to center them
    #Eagle Hill Rd high ground entrance
    gauges.append([1, -124.029176, 46.728031, 0, 1e+09])

    #Eagle Hill Rd and SR 105
    gauges.append([2, -124.027977, 46.726885, 0, 1e+09])

    #Tokeland Rd and SR 105
    gauges.append([3, -124.020344, 46.72461, 0, 1e+09])

    #Shoalwater Bay Gynasium
    gauges.append([4, -124.013834, 46.721922, 0, 1e+09])

    #Shoalwater Bay Gynasium Trail Middle
    gauges.append([5, -124.014507, 46.723688, 0, 1e+09])

    #Shoalwater Bay Gynasium Trail and SR 105
    gauges.append([6, -124.013470, 46.725850, 0, 1e+09])

    #Trail to Teal Duck Slough, closer to SR 105
    gauges.append([7, -124.010387, 46.727853, 0, 1e+09])

    #Trail to High Ground past Duck Slough Road and SR 105
    gauges.append([8, -124.004113, 46.727952, 0, 1e+09])

    #High Ground past Duck Slough Road (clearcut)
    gauges.append([9, -124.005254, 46.728384, 0, 1e+09])

    #Shoalwater Bay Tribal Court
    gauges.append([10, -124.016535, 46.721395, 0, 1e+09])

    #Shoalwater Bay Food Bank
    gauges.append([11, -124.0152, 46.721125, 0, 1e+09])

    #Tradewinds on the Bay
    gauges.append([12, -124.01304, 46.717688, 0, 1e+09])

    #Auntie Lee and Blackberry Lane entrance
    gauges.append([13, -123.990825, 46.710251, 0, 1e+09])

    #South Beach Regional fire Authority
    gauges.append([14, -123.995923, 46.711095, 0, 1e+09])

    #Shoalwater Museum
    gauges.append([15, -124.021153, 46.724253, 0, 1e+09])


if 0:
    lines = open('GHPC_Bridge_Coordinates.csv').readlines()
    gauges = []

    for line in lines[1:]:
        tokens = line.split(',')
        if (len(tokens) > 1):
            gaugeno=int(tokens[0])
            gauges.append([gaugeno, float(tokens[6]), float(tokens[5]), 0, 1e9])

gauges_array = np.array(gauges)

x_desired = gauges_array[:,1]
y_desired = gauges_array[:,2]

#FIX
#
onesixth = 1/(6*3600.)  # 1/6" used for shifting domain

x_edge = -136.75  - onesixth
y_edge =  38. - onesixth

# center in 1/3" cells:
dx = 1/(3.0*3600.)
dy = 1/(3.0*3600.)

x_centered, y_centered = center_points.adjust_xy(x_desired, y_desired,
                                                 x_edge, y_edge, dx, dy,
                                                 verbose=True)
gauges_array[:,1] = x_centered
gauges_array[:,2] = y_centered

for k in range(len(gauges)):
    # 3" flagregion.spatial_region = [-124.295,-123.8, 46.8,47.1]
    print('gauges.append([%i, %.8f, %.8f, %g, %g])' % tuple(gauges_array[k,:]))
