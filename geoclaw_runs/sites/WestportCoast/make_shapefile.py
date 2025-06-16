from pylab import *
import os

import fiona
import numpy as np
import rasterio
import rasterio.features
import rioxarray  # noqa
import shapely
import xarray as xr
from shapely import geometry
from shapely.geometry import Polygon, MultiPolygon, shape
from shapely.plotting import plot_polygon
from matplotlib import cm


def ds2shp(ds, fileout, varname, contours, epsg, geom_type="MultiLineString"):

    """
    Convert a xarray.Dataset to a shapefile, e.g. to make a shapefile
    for an inundation region that can be imported in GIS tools for plotting
    on a map.

    Written by Katy Barnhart

    Inputs:
        ds : xarray

        fileout : str
            Path of output vector file. If fileout ends with 'shp', a shapefile
            will be written. If it ends with 'geojson' a GeoJson will be written.
        varname : str
            Variable name to write out.
        contours : list
            A list of levels to use as the contour values.
        epsg : int
            EPSG code used to assign coordinate reference (optional)
        geom_type : str
            Either "MultiPolygon" for polygons or "MultiLineString" for line
            strings. The former gives the entire area whereas the latter gives
            the edges of the area impacted with varname>contours.
    Outputs:
        None
    """
    # create shapes
    shape_list = []
    for contour in contours:
        mask = ds[varname].values > contour

        # use rasterio.features to generate a polyline from a mask at each time.
        # first, get the polygon.
        shapes = [
            geometry.shape(feat)
            for feat, val in rasterio.features.shapes(
                mask.astype(np.int16), mask=mask, transform=ds.rio.transform(True)
            )
        ]

        # then, convert to polylines. get interiors and exteriors.
        if geom_type == "MultiLineString":
            out_rings = []
            for shape in shapes:
                out_rings.append(shape.exterior)
                out_rings.extend(shape.interiors)

            shape_list.append(
                {
                    "poly": shapely.MultiLineString(out_rings),
                    "props": {varname: float(contour)},
                }
            )
        elif geom_type == "MultiPolygon":
            shape_list.append(
                {
                    "poly": shapely.MultiPolygon(shapes),
                    "props": {varname: float(contour)},
                }
            )
        else:
            raise ValueError(
                f"digger.utils.fgtools.fg2shp given an invalid geom_type: {geom_type}"
            )
    # define ESRI schema, write each polygon to the file
    if "shp" in fileout:
        writer = "ESRI Shapefile"
    elif "geojson" in fileout:
        writer = "GeoJSON"
    else:
        raise ValueError("fileout must end with shp or geojson")

    schema = {
        "geometry": geom_type,
        "properties": {varname: "float"},
        "crs": f"EPSG:{epsg}",
    }
    with fiona.collection(fileout, "w", writer, schema) as output:
        for shape_info in shape_list:
            output.write(
                {
                    "properties": shape_info["props"],
                    "geometry": geometry.mapping(shape_info["poly"]),
                }
            )


def make_sample():
    """
    Define an array of water depths and then try to create a shapefile
    for the polygon where water depth > dry_tol
    Read this back in and plot it.
    """

    # make sample water depth data:
    x = linspace(-10,10,101)
    y = linspace(-10,10,51)[::-1]  # reverse order
    X,Y = meshgrid(x,y,indexing='xy')

    if True:
        # gaussians:
        h = exp(-0.1*(X**2 + Y**2)) + exp(-0.2*((X-4)**2 + (Y-4)**2)) \
            - exp(-0.8*((X-2.5)**2 + (Y-2.5)**2))  # with a hole
        h = maximum(h, 0) # so non-negative

    else:
        # use rectangles to get simpler polygons:
        h = where(logical_and(minimum(X,Y)>-4, maximum(X,Y)<6),1,0)
        h = where(logical_and(minimum(X,Y)>3, maximum(X,Y)<5),0,h)


    # convert to xarray dataset and plot:
    harray = xr.DataArray(h, dims=['y','x'], coords = {'x':x,'y':y})
    ds = xr.Dataset({'h':harray})
    ds.rio.set_spatial_dims(
                x_dim="x",
                y_dim="y",
                inplace=True,
            )
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    west=np.min(x)- dx/2
    east=np.max(x)+ dx/2
    south=np.min(y)+ dy/2
    north=np.max(y)- dy/2

    figure(1,figsize=(7,7))
    clf()
    clines = arange(0.1, 2, 0.1)
    imshow(harray,cmap=cm.Reds, extent=(west, east, south, north), interpolation='none')
    dry_tol = 0.1
    contour(X,Y,harray,levels=[dry_tol],colors='b',linewidths=1)
    axis('equal')

    # make a shapefile:
    outdir = 'shapefiles'
    os.system('mkdir -p %s' % outdir)
    fileout = '%s/sample.shp' % outdir
    varname = 'h'
    contours = [dry_tol]
    epsg = None
    ds2shp(ds, fileout, varname, contours, epsg,
                          geom_type="MultiPolygon")

    # read shapefile back in and plot:

    mp = MultiPolygon([shape(pol['geometry']) for pol in fiona.open(fileout)])

    if 0:
        mp = mp.simplify(0.4) # produces polygons with fewer vertices

    plot_polygon(mp, color='m', add_points=False)

    # print('Contents of MultiPolygon:')
    # for p in mp.geoms:
    #     print(p)

    savefig('compare.png')

if __name__ == '__main__':
    make_sample()
