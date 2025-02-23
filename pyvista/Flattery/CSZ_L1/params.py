"""
Parameters specific to the earthquake source (and the location).
"""

import numpy as np
import sys, os
import csv

from clawpack.geoclaw import fgout_tools, fgmax_tools
from clawpack.amrclaw.data import FlagRegion

#root_dir = os.path.abspath('../../..')  # should point to main BellinghamBayMaritime

#root_dir_Snoho2 = os.environ['Snoho2_SHARED']
#root_dir_Snoho2 = os.environ['Snoho2']
#root_dir = os.environ['WA_EMD_2020_SHARED']
root_dir = os.environ['WA_EMD_2020']
input_files_dir = root_dir + '/topo/input_files/'

#topodir_Snoho2 = os.path.join(root_dir_Snoho2, 'topo/topofiles')
#topodir = os.path.join(root_dir, 'topo/topofiles')
topodir = '/Users/rjl/topo/topofiles'
#print('Will look for some topo files in %s'  % topodir_Snoho2)
print('Will also look for some topo files in %s'  % topodir)
if not os.path.isdir(topodir):
    raise ValueError('Is this correct?  root_dir = %s' % root_dir)



#import data_FlagRegions
print('+++ fgmax_tools from ',fgmax_tools.__file__)

# -----------------------'------------------------------------
# Set loc to name of this fgmax point region, event to tsunami source:
loc = 'Flattery'
event = 'CSZ_L1'

print('Assuming loc = %s, event = %s' % (loc,event))

#Flattery
fgmax_extent = [-124.8, -124.52, 48.192, 48.435]
#Ozette
#fgmax_extent = [-124.8, -124.59, 48.0, 48.2]
#LaPush
#fgmax_extent = [-124.8, -124.4065, 47.765, 48.008]


restart = False    # if True, set proper checkpoint file in next line:
#restart_file = os.path.abspath('_output/fort.chkbbbbb')
restart_file = '/rc_scratch/rale6846/WA_EMD_2020/Runs/Flattery/CSZ_L1_v570/_output/fort.chkbbbbb'
# note needs to point to scratch directory on CU
# see the OUTPUT.txt file to get path to scratch directory


# ------------------
# setrun parameters:
# ------------------

output_style = 1
tfinal = 1.5*3600.
num_output_times = 1

t_start = 0*3600                   #for finest region and fgmax

arcsec16 = 1./(6*3600.)
arcsec13 = 1./(3*3600.)

lower = [0,0]
upper = [0,0]
num_cells = [0,0]

# Set the computational domain
lower[0] = -135.16 - arcsec16    # west longitude
upper[0] = -122.16- arcsec16    # east longitude
lower[1] = 38.5 - arcsec16       # south latitude
upper[1] = 53.5 - arcsec16       # north latitude

# Number of grid cells: Coarsest grid - 6' by 6' cells
num_cells[0] = 13
num_cells[1] = 15

dimensional_split = 'unsplit'  # 'unsplit' or 1 for dimensional_split

# AMR parameters:
# max number of refinement levels:
amr_levels_max = 7

# List of refinement ratios at each level (length at least mxnest-1)
# dx = dy = 1deg, 6', 2', 30", 6", 2", 1/3":
refinement_ratios = [10,3,4,5,3,6]

topofiles = []

# 1-minute topo:
topofiles.append([3, topodir + '/etopo1_-163_-122_38_63.asc'])

# 2-second topo:
topofiles.append([3,  topodir + '/PT_2sec_center.asc'])
topofiles.append([3,  topodir + '/SJdF_2sec_center.asc'])
topofiles.append([3,  topodir + '/NWOuterCoast_2s_mhw_lapush.asc'])

# 1/3-second topo needed around fgmax region for all sources:
topofiles.append([3,  topodir + '/NWOuterCoast_13_mhw_merged.asc'])

# dtopo file
dtopo_path =   '/Users/rjl/B/dtopo/dtopofiles/CSZ_L1-extended-pmel.tt3'
dtopofiles = [[3,dtopo_path]]

# regions:
regions = []

#tstart_regions = 0*60.

# Levels 1 to 7, dx = dy = 1deg, 6', 2', 30", 6", 2", 1/3":

#whole domain, require 1deg, but allow 2'
#Just keeping the second 1,3 domain after 3 hours from run completion.

regions.append([1,3,0.,3*3600.,lower[0]-0.1,upper[0]+0.1,lower[1]-0.1,upper[1]+0.1])
regions.append([1,1,3*3600.,1e9,lower[0]-0.1,upper[0]+0.1,lower[1]-0.1,upper[1]+0.1])
regions.append([1,3,3*3600.,1e9,-131.0,-122.06,45.89,51.10])

#source region
regions.append([3,4,0.,1*30.,-135.0,-122.0,39.0,53.0])

# Close to the entrance to the Strait along the coast, require 2 min ocean, but allow 30"
# in the Ruled Rectangle called Coast_46_51 below, turned on at time 0.

#This is the mouth of the Strait up to Port Angeles, roughly, require 2', but allow 6".
#Followed by an extension of this south to include the LaPush area
#Turn on at 0 hours
regions.append([3,5,0*3600.,1e9,-125.115,-123.56+0.005,48.045,48.715])
regions.append([3,5,0*3600.,1e9,-125.05,-124.35,47.55,48.04])

#These are rectangles around the study area, require 30", could be 2"
regions.append([4,6,0*3600.,1e9,-124.888,-124.51,47.9501,48.68])
regions.append([4,6,0*3600.,1e9,-124.888,-124.436,47.75,47.9501])

#This is PA eastward in the Strait, getting closer to Whidbey Island,
#so do 3,5, require 2', but allow 6". If further north above 48.395, dont refine as much now.

# Turning on sooner to get the negative wave
regions.append([3,5,0.0*3600.,1e9,-123.56+.005,-122.91,48.0,48.395])
regions.append([3,4,0.0*3600.,1e9,-123.56+.005,-122.91,48.395,48.785])
regions.append([3,5,10.*60.,1e9,-122.91,-122.66,48.0,48.395])
regions.append([3,4,10.*60.,1e9,-122.91,-122.42,48.395,48.785])

#This is northward around Orcas and northern water, first is 6" allowed,
#    30" required, and the further north one is 30" required, 6" allowed also.
#regions.append([3,4,tstart_regions,1e9,-124.01,-122.482,48.833,49.585])
#Dont need these for the Outer coast project
#regions.append([3,4,tstart_regions,1e9,-124.01,-122.482,48.785,49.585])
#regions.append([3,4,tstart_regions,1e9,-125.24,-124.01,49.23,50.125])

#######We used 47.851 for Bellingham Maritime CSZ, but our south study region close to here
#######   so give a slightly bigger 2sec region to 47.914.
#Changed the 47.851 in the 3,6 region below to now be 47.914
#x2_thin = -122.50 + .2
#regions.append([4,6,tstart_regions,1e9,-122.56,x2_thin,47.73,47.914])
#regions.append([3,5,t_start,1e9,-122.56,x2_thin,47.01,47.73])
##regions.append([3,5,tstart_regions,1e9,-122.56,x2_thin,47.01,47.914])

# Flag regions defined as RuledRectangles:
RRdir = os.path.join(root_dir, 'topo/RuledRectangles')

flagregions = []

# append as many flagregions as desired to this list:

flagregion = FlagRegion(num_dim=2)
flagregion.name = 'Region_Coast_46_51'
flagregion.minlevel = 3
flagregion.maxlevel = 4
flagregion.t1 = 0.*3600.
flagregion.t2 = 1e9
flagregion.spatial_region_type = 2  # Ruled Rectangle
flagregion.spatial_region_file = os.path.abspath(RRdir + \
        '/RuledRectangle_Coast_46_51.data')
flagregions.append(flagregion)

flagregion = FlagRegion(num_dim=2)
flagregion.name = 'Flattery'
flagregion.minlevel = 7
flagregion.maxlevel = 7
flagregion.t1 = t_start
flagregion.t2 = 1e9
flagregion.spatial_region_type = 2  # Ruled Rectangle
flagregion.spatial_region_file = os.path.abspath(input_files_dir + \
        '/RuledRectangle_Flattery.data')
flagregions.append(flagregion)

flagregion = FlagRegion(num_dim=2)
flagregion.name = 'Ozette'
flagregion.minlevel = 5
flagregion.maxlevel = 6
flagregion.t1 = t_start
flagregion.t2 = 1e9
flagregion.spatial_region_type = 2  # Ruled Rectangle
flagregion.spatial_region_file = os.path.abspath(input_files_dir + \
        '/RuledRectangle_Ozette.data')
flagregions.append(flagregion)

flagregion = FlagRegion(num_dim=2)
flagregion.name = 'LaPushN'
flagregion.minlevel = 5
flagregion.maxlevel = 6
flagregion.t1 = t_start
flagregion.t2 = 1e9
flagregion.spatial_region_type = 2  # Ruled Rectangle
flagregion.spatial_region_file = os.path.abspath(input_files_dir + \
        '/RuledRectangle_LaPush_N.data')
flagregions.append(flagregion)

flagregion = FlagRegion(num_dim=2)
flagregion.name = 'LaPushS'
flagregion.minlevel = 5
flagregion.maxlevel = 6
flagregion.t1 = t_start
flagregion.t2 = 1e9
flagregion.spatial_region_type = 2  # Ruled Rectangle
flagregion.spatial_region_file = os.path.abspath(input_files_dir + \
        '/RuledRectangle_LaPush_S.data')
flagregions.append(flagregion)

# Tolerance for refining if not otherwise constained:
wave_tolerance = 0.05

#--------------
# fgmax points:

# NOTE new format, now specify fgmax_grids

fg = fgmax_tools.FGmaxGrid()
fg.point_style = 4       # specified using mask in topo_type 3 file

if 1: #none yet
    fgmax_grids = []
if 0:
    fg.xy_fname = input_files_dir + 'fgmax_pts_Flattery.data'
    fg.min_level_check = amr_levels_max  # which levels to monitor max on
    fg.tstart_max = t_start + 120
    fg.tend_max = 1.e10              # when to stop monitoring max values
    fg.dt_check = 3.                # how often to update max values
    fgmax_grids = [fg]

#--------------
# fgout grid:

fgout_grids = []

fgout_extent = [-124.76, -124.55, 48.29, 48.41]
dx_fgout = 1 * 1/3600.  # degrees
dy_fgout = dx_fgout
dt_fgout = 30  # seconds

fgout = fgout_tools.FGoutGrid()
fgout.fgno = 3
fgout.point_style = 2  # uniform rectangular x-y grid
fgout.output_format = 'binary32'
fgout.x1 = fgout_extent[0]
fgout.x2 = fgout_extent[1]
fgout.y1 = fgout_extent[2]
fgout.y2 = fgout_extent[3]
fgout.nx = int(round((fgout.x2 - fgout.x1)/dx_fgout))
fgout.ny = int(round((fgout.y2 - fgout.y1)/dy_fgout))
fgout.tstart = 0.
#fgout.tstart = 0.  # SHORT TEST
fgout.tend = tfinal
fgout.nout = int(round(((fgout.tend - fgout.tstart)/dt_fgout))) + 1
#q_out_vars = [1,2,3,4] # h,hu,hv,eta for shallow code
q_out_vars = [1,4] # h,eta for shallow code
fgout.q_out_vars = q_out_vars
fgout_grids.append(fgout)


#--------------
# gauges:

# Add any gauges from the file provided by DNR of desired locations
# for gauges that fall in the region being modeled.
#
gauges = []

#tstart_gauges = 0.                  # when to start collecting gauges
tstart_gauges =  t_start  + 120      # start 120sec after finest grid on

# list of possible gauges:
gauges_csv_file = root_dir + '/info/tide_gauge_locations.csv'

gaugenos = []   # gauges to always include

f = open(gauges_csv_file,'r')
csv_reader = csv.reader(f)
gaugeloc = {}
gaugex = {}
gaugey = {}
for row in csv_reader:
    try:
        xg = float(row[2])
        yg = float(row[3])
    except:
        continue  # go to next row
    #gaugeno = int(row[0]) + 1000
    gaugeno = int(row[0])
    gaugeloc[gaugeno] = row[4]
    gaugex[gaugeno] = xg
    gaugey[gaugeno] = yg

    if ((fgmax_extent[0] <= xg <= fgmax_extent[1]) and \
        (fgmax_extent[2] <= yg <= fgmax_extent[3])) and \
        gaugeno not in gaugenos:
            print('+++ adding gauge ',gaugeno)
            gaugenos.append(gaugeno)  # include this gauge

f.close()

for gaugeno in gaugenos:
    gauges.append([gaugeno, gaugex[gaugeno], gaugey[gaugeno], tstart_gauges, 1e9])
    print('Including gauge %3i at %s' % (gaugeno, gaugeloc[gaugeno]))


# gauges is now a list of lists, possibly empty
# append other gauges if desired, via e.g.
#  gauges.append([gaugeno, x, y, t1, t2])

#--------------

# subsidence in study region only for CSZ_L1 or CSZ_XL1 event:
#y1_fgmax = fgmax_extent[2]

variable_eta_init = (('CSZ' in event) or (loc == 'Ozette'))
print('+++ variable_eta_init = ',variable_eta_init)

topo_missing = 100.   # shouldn't be any missing topo

#--------------
# force_dry_init:
# If there are regions below MHW that should be intially dry,
# the force_dry_init file indicates such points.

force_dry_file_used = False
force_dry_fname = input_files_dir + '/force_dry_init.data'

# time period over which these points should be forced dry when new finest
# level grids created. Should be slightly after grids are first created,
# and before tsunami starts inundating:
force_dry_tend = t_start + 120.



# ------------------
# setplot parameters:
# ------------------


plot_shore = False   # still need to generate file with shorelines!
if plot_shore:
    shorelines = np.loadtxt(topodir + '/SkagitIsland_shoreline.txt')
else:
    shorelines = None

zoom_extent = fgmax_extent
figsize = (6,8) # for zoom extent
