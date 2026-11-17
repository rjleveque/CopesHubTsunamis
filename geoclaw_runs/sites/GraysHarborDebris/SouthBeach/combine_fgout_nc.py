"""
Read in the fgout output from a GeoClaw run and combine all frames
into a single netCDF file.

This is a simple wrapper for fgout_tools.write_netcdf

Set to capture the following quantities of interest: h, u, v, B
"""

from clawpack.geoclaw import fgout_tools

loc = 'SouthBeach'
event = 'BL10S'
#event = 'BR10M'
outdir = f'_output_{event}'

fgno = 1
output_format = 'binary32'
fname_nc=f'{loc}_fgout_arrays_{event}.nc'

fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format)
fgout_grid.read_fgout_grids_data()

fgout_frames = []
for fgframeno in range(1,fgout_grid.nout+1):
    fgout = fgout_grid.read_frame(fgframeno)
    fgout_frames.append(fgout)
    
fgout_tools.write_netcdf(fgout_frames, fname_nc=fname_nc,
                         qois = ['h','u','v','B'], datatype='f4',
                         include_B0=False, include_Bfinal=False,
                         description='', verbose=True)

