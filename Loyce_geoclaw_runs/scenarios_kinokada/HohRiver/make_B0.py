from pylab import *
import os,sys
from clawpack.geoclaw import topotools
from clawpack.geoclaw import fgmax_tools

### There are 3 fgmax files for the Copes Hub Runs
#fgname=['Hoh_B0.asc','Ruby_B0.asc','Kalaloch_B0.asc']

### Made a new Hoh file for L1, as it needed to be bigger
fgname=['Hoh_forL1_B0.asc']

outdir0 = '_output_B0_forL1'

print('Using output from outdir0 = ',outdir0)

# Read fgmax data:
fgmax_input_file_name = outdir0 + '/fgmax_grids.data'
print('fgmax input file: \n  %s' % fgmax_input_file_name)

###### Loop over the fgmax grids
#for fgno in [1,2,3]:
for fgno in [1]:
    fg0 = fgmax_tools.FGmaxGrid()
    fg0.read_fgmax_grids_data(fgno=fgno, data_file=fgmax_input_file_name)

    # read fgmax results with indexing 'xy' so B0[j,i] is for x[i],y[j]
    # and so B0 has the orientation expected for a topotools.Topography
    fg0.read_output(outdir=outdir0, indexing='xy')

    B0 = where(fg0.B > -1e9, fg0.B, -9999.)

    print('retrieved B0 using fgmax_tools from clawpack of shape ',B0.shape)
    print('shape of fg0.x was: ',shape(fg0.x))
    print('shape of fg0.y was: ',shape(fg0.y))
    print('plan to execute: topoB0.set_xyZ(fg0.x, fg0.y, B0)')
    print(' ')


    fname = 'input_files/' + fgname[fgno-1]
    topoB0 = topotools.Topography()
    topoB0.set_xyZ(fg0.x, fg0.y, B0)
    topoB0.write(fname, topo_type=3, header_style='asc', Z_format='%.3f')

    print('shape of fg0.x was: ',shape(fg0.x))
    print('shape of fg0.y was: ',shape(fg0.y))
    print('executed topoB0.set_xyZ(fg0.x, fg0.y, B0): ')
    print('saved %s using topotools with B0 of shape ' % fname, B0.shape)
    print('shape of topoB0.Z was: ',shape(topoB0.Z))
    print('which hopefully matches shape of B0')
    print(' ')

