"""
Try out fg2shp function
"""

from pylab import *
from fg2output import _get_ds
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import colormaps
from clawpack.visclaw import plottools
import make_shapefile
import geopandas
from shapely.plotting import plot_polygon
import shapely
import fiona


outdir = 'CSZ_M1/_output'
fgno = 1

filein = '%s/fgmax%s.txt' % (outdir, str(fgno).zfill(4))
fileout = 'CSZ_M1_inundation.shp'
varname = 'h_max'
contours = [0.02]
epsg = 32610
fgtype = 'fgmax'
geom_type: str = "MultiLineString"

ds = _get_ds(filein, fgtype, epsg)
make_shapefile.ds2shp(ds, fileout, varname, contours, epsg,
                      geom_type="MultiPolygon")

#ds = fg2output.fg2shp(filein,fileout,varname,contours,epsg,fgtype,
#                      geom_type, maskB0)

# now try reading it back in and plotting:
if 0:
    inundation_region = geopandas.read_file(fileout)
else:
    polys = [shapely.geometry.shape(pol['geometry']) \
                                     for pol in fiona.open(fileout)]
    inundation_region = shapely.geometry.MultiPolygon(polys[0])


fg = fgmax_tools.FGmaxGrid()
fgmax_input_file_name = outdir + '/fgmax_grids.data'
print('fgmax input file: \n  %s' % fgmax_input_file_name)
fg.read_fgmax_grids_data(fgno=fgno, data_file=fgmax_input_file_name)

fg.read_output(outdir=outdir)

fg.B0 = fg.B

zmin = -10.
zmax = 10.
land_cmap = colormaps.make_colormap({ 0.0:[0.1,0.4,0.0],
                                     0.25:[0.0,1.0,0.0],
                                      0.5:[0.8,1.0,0.5],
                                      1.0:[0.8,0.5,0.2]})

sea_cmap = colormaps.make_colormap({ 0.0:[0,0,1], 1.:[.8,.8,1]})

cmap, norm = colormaps.add_colormaps((land_cmap, sea_cmap),
                                     data_limits=(zmin,zmax),
                                     data_break=0.)

fig,axs = subplots(1,2,figsize=(12,8))
ax = axs[0]
#plottools.pcolorcells(fg.X, fg.Y, fg.B, ax=ax, cmap=cmap, norm=norm)
inundated = ma.masked_where(logical_or(fg.B<0, fg.h<0.01), fg.h)
pc = ax.pcolormesh(fg.X,fg.Y,inundated)
colorbar(pc, shrink=0.7)
ax.contour(fg.X,fg.Y,fg.B,[0],colors='g', linewidths=0.9)
ax.set_aspect(1/cos(pi*fg.Y.mean()/180))
ax.set_xlim(-124.16,-124.06)
ax.set_ylim(46.75,46.92)
xticks(rotation=20);

ax = axs[1]
plot_polygon(inundation_region, ax=ax, color='m', add_points=False)
ax.contour(fg.X,fg.Y,fg.B,[0],colors='g', linewidths=0.9)
ax.set_aspect(1/cos(pi*fg.Y.mean()/180))
ax.set_xlim(-124.16,-124.06)
ax.set_ylim(46.75,46.92)
xticks(rotation=20);
