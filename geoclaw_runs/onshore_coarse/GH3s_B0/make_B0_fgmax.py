
from numpy import save
from clawpack.geoclaw import fgmax_tools, topotools

location = 'GH3s'
fgno = 1

fgmax = fgmax_tools.FGmaxGrid()
fgmax.read_fgmax_grids_data(fgno)
fgmax.read_output(indexing='xy')   # read with topofile orientation [j,i]

B0 = fgmax.B
if B0.mask.sum() > 0:
    print('*** Warning B0 has some masked values')

topoB0 = topotools.Topography()
topoB0.set_xyZ(fgmax.x, fgmax.y, B0)

fname = f'{location}_fgmax{fgno}_B0.asc'
topoB0.write(fname, topo_type=3, header_style='asc', Z_format="%.3f")
print('Created ',fname)

