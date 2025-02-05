
from pylab import *
import os
from clawpack.geoclaw import fgmax_tools
from clawpack.geoclaw import topotools

sys.path.insert(0,'../common_code')
#import fgmax_tools  # uses local version with transposed arrays
                    # should appear in v5.10.0

outdir0 = os.path.abspath('Newport_B0/_output')
print('Using output from outdir0 = ',outdir0)

# Read fgmax data:
fg0 = fgmax_tools.FGmaxGrid()
fgmax_input_file_name = outdir0 + '/fgmax_grids.data'
print('fgmax input file: \n  %s' % fgmax_input_file_name)
fg0.read_fgmax_grids_data(fgno=1, data_file=fgmax_input_file_name)
fg0.read_output(outdir=outdir0, indexing='xy')  # xy to conform to topofile
B0 = where(fg0.B > -1e9, fg0.B, -9999.)

fname = 'fgmax0001_13s_B0.asc'
topoB0 = topotools.Topography()
topoB0.set_xyZ(fg0.x, fg0.y, B0)
topoB0.write(fname, topo_type=3, header_style='asc', Z_format='%.3f')
print('saved %s with B0 of shape ' % fname, B0.shape)

# Read fgmax data:
fg0 = fgmax_tools.FGmaxGrid()
fgmax_input_file_name = outdir0 + '/fgmax_grids.data'
print('fgmax input file: \n  %s' % fgmax_input_file_name)
fg0.read_fgmax_grids_data(fgno=2, data_file=fgmax_input_file_name)
fg0.read_output(outdir=outdir0, indexing='xy')  # xy to conform to topofile
B0 = where(fg0.B > -1e9, fg0.B, -9999.)

fname = 'fgmax0002_13s_B0.asc'
topoB0 = topotools.Topography()
topoB0.set_xyZ(fg0.x, fg0.y, B0)
topoB0.write(fname, topo_type=3, header_style='asc', Z_format='%.3f')
print('saved %s with B0 of shape ' % fname, B0.shape)
