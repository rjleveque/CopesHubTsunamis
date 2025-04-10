from pylab import *
import pyvista as pv
from clawpack.geoclaw import topotools
from clawpack.clawutil.data import get_remote_file
import os


global etamesh

topodir = '/Users/rjl/topo/topofiles/'
topo_fname = topodir + 'Flattery_13s_mhw.asc'
if not os.path.isfile(topo_fname):
    # need to download it first...
    url = 'http://depts.washington.edu/clawpack/geoclaw/topo/WA/' + topo_fname
    get_remote_file(url, output_dir='.', file_name=topo_fname, verbose=True)

topo = topotools.Topography(topo_fname)

extent = [-124.76, -124.55, 48.29, 48.41]
topo = topo.crop(extent)

# matplotlib topo plot:
#topo.plot(limits=(-25,25))

# background image to use a texture on warped topo:
bg_file = 'Flattery_img.jpg'

if not os.path.isfile(bg_file):
    print('First run Quillayute_bg_image.py or provide other bg image file')
    print('Settting bg_file to None')
    bg_file = None

bg_file = None  # for no background image, topo colored by elevation instead

z = array([0.])
x = (topo.x - topo.x[0]) * 111e3 * cos(topo.y.mean()*pi/180)
y = -(topo.y - topo.y[0]) * 111e3
print('xmax = %.1fm, ymax = %.1fm' % (x.max(),y.max()))
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

B = flipud(topo.Z)
Bmax = 500. #150.
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='C')

warpfactor = 3   # amplification of elevations

topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

make_snapshots = False  # True to capture a few frames, False for slider bar

p = pv.Plotter(off_screen=make_snapshots)

if bg_file:
    # Add bg image as texture:
    bg_extent = extent
    texture = pv.read_texture(bg_file)

    x1,x2 = (asarray(bg_extent[:2]) - topo.x[0]) * 111e3 * cos(topo.y.mean()*pi/180)
    y1,y2 = (asarray(bg_extent[2:]) - topo.y[0]) * 111e3

    o = (x1, y1, 0.)
    u = (x2, y1, 0.)
    v = (x1, y2, 0.)

    mapped_surf = topowarp.texture_map_to_plane(o, u, v)
    p.add_mesh(mapped_surf,texture=texture)

else:
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-0.8*Bmax, 0.8*Bmax))


sea_level = 0.
eta = where(B < sea_level, sea_level, nan)
topoxyz.point_data['eta'] = eta.flatten(order='C')
etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
etamesh = p.add_mesh(etawarp,color='c')


p.camera_position =   [(7731.396, -15682.454, 3888.938), (7806.941, -3573.918, -2200.799), (0.006, 0.449, 0.893)]

def set_sea_level(sea_level):
    global etamesh
    eta = where(B < sea_level, sea_level, nan)
    topoxyz.point_data['eta'] = eta.flatten(order='C')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)
    etamesh = p.add_mesh(etawarp,color='c')

    if not make_snapshots:
        # round off entries in p.camera_position and print out, so user
        # can copy and paste into this script once good position is found:
        camera_position = list(p.camera_position)
        for i,a in enumerate(camera_position):
            b = []
            for j in range(len(a)):
                b.append(round(a[j],3))
            camera_position[i] = tuple(b)
        print('p.camera_position = ', camera_position)

if not make_snapshots:
    p.add_title('MHW after sea level rise / subsidence')
    p.add_slider_widget(set_sea_level, [-5,5], value=0, title='Sea Level',
                        pointa=(0.1,0.1), pointb=(0.4,0.1),)
    p.show(window_size=(2500,1500))
else:
    p.window_size = (2500,1500)
    for slr in [0,1,2,3]:
        set_sea_level(slr)
        p.add_title('MHW after %i m subsidence (or sea level rise)' % slr)
        fname_png = 'Flattery_mhw%im.png' % slr
        p.screenshot(fname_png)
        print('Created ',fname_png)
    p.close()
