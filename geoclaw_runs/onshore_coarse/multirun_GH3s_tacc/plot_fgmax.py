
"""
Simple demo to load fgmax results for one event, plot inundation depth onshore.
"""

from pylab import *
import os
from numpy import ma # masked arrays
from matplotlib import colors
from clawpack.visclaw import plottools
from clawpack.clawutil.data import ClawData

def load_fgmax(outdir, fgno, fname_B0=None):
    from clawpack.geoclaw import fgmax_tools, topotools
    # Read fgmax data:
    fg = fgmax_tools.FGmaxGrid()
    fgmax_input_file_name = outdir + '/fgmax_grids.data'
    print('fgmax input file: \n  %s' % fgmax_input_file_name)
    fg.read_fgmax_grids_data(fgno=fgno, data_file=fgmax_input_file_name)

    # determine time interval used for fgmax:
    clawdata = ClawData()
    clawdata.read(os.path.join(outdir, 'claw.data'), force=True)
    try:
        t_hours = min(clawdata.tfinal, fg.tend_max) / 3600.
    except:
        t_hours = nan

    fg.read_output(outdir=outdir, indexing='xy')  # so array layout same as topofile

    #### Read pre-seismic B0 from special run with no dtopo specified

    if fname_B0:
        topoB0 = topotools.Topography()
        topoB0.read(fname_B0, topo_type=3)
        B0 = topoB0.Z
        B0_masked = ma.masked_array(B0, fg.B.mask)
        fg.B0 = B0_masked
    else:
        fg.B0 = fg.B
        print('Assuming B0 = fg.B')

    return fg, t_hours


# load fgmax results from particular event:

outdir = 'geoclaw_outputs_tacc/_output_BL10D'  # or full path to output directory
fgno = 1  # fixed grid number

# set path to B0 topography on fgmax grid computed via run with no dtopo:
fname_B0 = '../GH3s_B0/GH3s_fgmax1_B0.asc'

# or set to None  to use fgmax.B to determine shoreline and what's onshore
# (but fgmax.B captured during or after the rupture, may include dz)
#fname_B0 = None

fgmax, t_hours = load_fgmax(outdir, fgno, fname_B0)
print(f'loaded fgmax results over simulated time {t_hours} hours')

# Use B0 to determine "onshore" points:
h_onshore = ma.masked_where(fgmax.B0 < 0, fgmax.h)
h_wet_onshore = ma.masked_where(h_onshore==0., h_onshore)


# simple plot for illustraton:
fig,ax = subplots(figsize=(8,8))

# colormap for depth:
bounds_depth = array([1e-6,1,2,4,6,10,12])

cmap_depth = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],
                 [1,.7,.7], [1,.4,.4], [1,0,0]])

# Set color for value exceeding top of range to purple:
cmap_depth.set_over(color=[1,0,1])

norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

pc = plottools.pcolorcells(fgmax.X, fgmax.Y, h_wet_onshore,
                           cmap=cmap_depth, norm=norm_depth)

colorbar(pc, extend='max', shrink=0.7, label='meters')

# pre-seismic shoreline:
contour(fgmax.X, fgmax.Y, fgmax.B0, [0], colors='g', linewidths=0.7)

ax.set_aspect(1/cos(pi*47/180))
