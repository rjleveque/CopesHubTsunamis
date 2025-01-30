"""
Using tiles downloaded from 
https://www.ngdc.noaa.gov/thredds/catalog/tiles/nthmp/tiled_19as/catalog.html

(Reading directly from url's did not seem to work)

"""

from clawpack.geoclaw import topotools

topodir = '/Users/rjl/topo/CUDEM'

extent = [-124.15, -124, 44.55, 44.68]
url = 'https://www.ngdc.noaa.gov/thredds/dodsC/tiles/nthmp/tiled_19as/ncei19_n44x75_w124x25_navd88_2021.nc'
fname = 'ncei19_n44x75_w124x25_navd88_2021.nc'
topo1 = topotools.read_netcdf(os.path.join(topodir,fname), extent=extent,
                             coarsen=1, verbose=True)


extent = [-124, -123.9, 44.55, 44.68]
url = 'https://www.ngdc.noaa.gov/thredds/dodsC/tiles/nthmp/tiled_19as/ncei19_n44x75_w124x00_navd88_2021.nc'
fname = 'ncei19_n44x75_w124x00_navd88_2021.nc'
topo2 = topotools.read_netcdf(os.path.join(topodir,fname), extent=extent,
                             coarsen=1, verbose=True)


# crop 1/9" around Newport and Yaquina Bay entrance:

extent = [-124.09, -124.0, 44.6, 44.66]
topo1a = topo1.crop(filter_region=extent)
fname = 'Newport_19s.asc'
topo1a.write(fname,topo_type=3,header_style='asc',Z_format='%12.5e')
print('Created ',fname)

# coarsen to 1/3"

topo1b = topo1.crop(coarsen=3)
fname = 'Newport_13s.asc'
topo1b.write(fname,topo_type=3,header_style='asc',Z_format='%12.5e')
print('Created ',fname)

topo2b = topo2.crop(coarsen=3)
fname = 'NewportE_13s.asc'
topo2b.write(fname,topo_type=3,header_style='asc',Z_format='%12.5e')
print('Created ',fname)

fig,ax = subplots()
topo1b.plot(limits=[-50,50],axes=ax,
        cb_kwargs={'extend':'both', 'shrink':0.7})
topo2b.plot(limits=[-50,50],axes=ax,add_colorbar=False)
xlim(-124.15, -123.9)
plot([-124,-124],[44.55,44.68],'k',linewidth=0.7)
title('Newport_13s and NewportE_13s')
fname = 'Newport_13s.png'
savefig(fname, bbox_inches='tight')
print('Created ',fname)

fig,ax = subplots()
topo1a.plot(limits=[-50,50],axes=ax,
        cb_kwargs={'extend':'both', 'shrink':0.7})
title('Newport_19s')
fname = 'Newport_19s.png'
savefig(fname, bbox_inches='tight')
print('Created ',fname)

if 0:
    # VDATUM didn't work on this file...
    # Create .txt files for VDATUM conversion:
    extent = [-124.07, -124.035, 44.615, 44.635]
    topo1c = topo1b.crop(filter_region=extent)
    fname = 'Newport_13s_marina.txt'
    topo1c.write(fname,topo_type=1,Z_format='%12.5e')
    print('Created ',fname)
    topo1c_mhw = topotools.Topography()
    topo1c_mhw.read('Newport_13s_marina_result.txt', topo_type=1)

    fig,ax = subplots()
    topo1c_mhw.plot(limits=[-50,50],axes=ax,
            cb_kwargs={'extend':'both', 'shrink':0.7})
    title('Newport_13s_marina')
    fname = 'Newport_13s_marina.png'
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)

# VDATUM at a few selected points shows MHW approx NAVD88 - 2.05 m
# https://vdatum.noaa.gov/vdatumweb/vdatumweb?a=191901220200922#Point


topos = [(topo1a,'Newport_19s_mhw'),
         (topo1b,'Newport_13s_mhw'),
         (topo2b,'NewportE_13s_mhw'),]

for (topo,name) in topos:
    Z_mhw = topo.Z - 2.05
    topo_mhw = topo.crop()
    topo_mhw.set_xyZ(topo.x, topo.y, Z_mhw)
    fig,ax = subplots()
    topo_mhw.plot(limits=[-50,50],axes=ax,
            cb_kwargs={'extend':'both', 'shrink':0.7})
    title(name)
    fname = '%s.asc' % name
    topo.write(fname,topo_type=3,header_style='asc',Z_format='%12.5e')
    print('Created ',fname)
