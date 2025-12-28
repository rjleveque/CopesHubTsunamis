"""
Testing the relevant asce gauges to see if the GeoClaw
output is compliant. Will use Westport since this will be
just pre to 2028 being adopted and we want to know some things!
The gauge locations changed from before.

 Westport 2025, Gauges 98 to 143
 Westport 2028, Gauges 127 to 164 

"""

# for making plots on remote machine:
import matplotlib
matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os
import copy
import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import dtopotools
import interptools
import params  #sets root_dir, loc, event, fgmax_extent

# SETTING directories where input data is found:
rundir = os.getcwd()
root_dir        = params.root_dir
assert os.path.isdir(root_dir), '*** did not find params.root_dir = %s' \
        % root_dir

out_dir         = params.out_dir
assert os.path.isdir(out_dir), '*** did not find output directory = %s' \
        % out_dir

plot_dir        = params.plot_dir
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

other_figures_dir    = plot_dir + '/_other_figures'
if not os.path.isdir(other_figures_dir):
    os.mkdir(other_figures_dir)

dtopo_dir = params.dtopo_dir
source_name = params.event
dtopofile = dtopo_dir + '/' + source_name + '.dtt3'
dtopotype = 3  # format of dtopo file

### 2028
asce_powell_dir = root_dir + '/info/asce_powell_gauges'

## Directories used
print('rundir = %s' % rundir)
print('root_dir is:  ',root_dir)
print('out_dir is: ',out_dir)
print('plot_dir is: ',plot_dir)
print('other_figures_dir is: ',other_figures_dir)
print('dtopo_dir is: ',dtopo_dir)
print('dtopo file is: ',dtopofile)
print('asce_powell_dir is: ',asce_powell_dir)
print(' ')


### 2028, Loyces file below has the lines from north (starting at 1) to south (ending at 551) 
###       from north to south increasing numbering.  ID[0]=1, ID[-1]=551.  The ID has the gaugeno
###       and it is always one more than the index number into the arrays ID, lon, lat, etc.

#pathname = asce_powell_dir + '/7-28_loyce_north-to-south-32km.txt'
pathname = '7-28_buried-locking_str10-deep.csv'
ID,lon,lat = np.loadtxt(pathname, delimiter=',',\
                                  skiprows=1,usecols=(0,1,2),unpack=True)
print(' ID[0],lon[0],lat[0] are:', ID[0], lon[0], lat[0])
print(' ')
