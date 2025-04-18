from pylab import *
import os
import pyvista as pv
from clawpack.geoclaw import fgout_tools
from clawpack.visclaw import animation_tools
from clawpack.visclaw import colormaps, geoplot
from time import sleep
from clawpack.geoclaw import topotools


rundir = os.getcwd()
event = os.path.split(rundir)[-1]

#event = 'buried-locking-mur13-deep'

if 0:
    topo = topotools.Topography('/Users/rjl/topo/topofiles/GH_tiles_2021_13s.asc')
    extent = [-124.185, -124.1325, 46.9725, 46.99]  # OS2 (fgout grid)
    topo = topo.crop(extent)

cmap_water = colormaps.make_colormap({0:'c', 1:[0,0.5,1]})

fgno = 1  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output

# where to save a png file for each frame, for constructing an animation:
framedir = '_frames'
os.system('mkdir -p %s' % framedir)
os.system('rm %s/*' % framedir)  # remove frames from previous version

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format)
fgout_grid.read_fgout_grids_data(fgno)
#fgout_grid.read_fgout_grids_data_pre511(fgno)
fgout = fgout_grid.read_frame(10) # after subsidence

# convert x,y to meters, roughly:

x_mean = fgout.x.mean()
y_mean = fgout.y.mean()
x = (fgout.x - x_mean) * 111e3 * cos(y_mean*pi/180)
y = -(fgout.y - y_mean) * 111e3
z = array([0.])
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

warpfactor = 1  # amplification of elevations
Bmax = 25.  # truncate hills at this elevation

def flipZ(Z):
    return fliplr(Z)

B0 = fgout.B  # save this as B0 for later plots of h+B0
B = flipZ(B0)
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='F')
topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

#import pdb; pdb.set_trace()

extent = [-124.14, -124.07, 46.88, 46.92]  # WestportN

GE_file = '../GE_WestportN.jpg'
GE_extent = extent


if 1:
    texture = pv.read_texture(GE_file)
    fname_mp4 = 'WestportN_animation.mp4'

#texture = pv.read_texture('./GE_Hoquiam_fgout2.jpg')
GE_extent = fgout.extent_edges  # must agree for texture mapping
x2,x1 = (asarray(GE_extent[:2]) - x_mean) * 111e3 * cos(y_mean*pi/180)
y2,y1 = (asarray(GE_extent[2:]) - y_mean) * 111e3

if 0:
    o = (x1, y1, 0.)
    u = (x2, y1, 0.)
    v = (x1, y2, 0.)

o = (x2, y2, 0.)
u = (x1, y2, 0.)
v = (x2, y1, 0.)

mapped_surf = topowarp.texture_map_to_plane(o, u, v)

# replace land surface by nan in eta so it only shows water:
eta = where(fgout.h>0.1, fgout.eta, nan)
eta = flipZ(eta)

make_animation = True

#p = pv.Plotter(off_screen=make_animation, lighting='none')
p = pv.Plotter(off_screen=make_animation, lighting='three lights')
#p = pv.Plotter(off_screen=make_animation)

p.window_size=(2300,1300)

if 0:
    light = pv.Light()
    light.set_direction_angle(30, 20)
    p.add_light(light)
    light = pv.Light()
    light.set_direction_angle(30, 0)
    p.add_light(light)

if 1:
    light1 = pv.Light(
        position=(-1000,0,200),
        focal_point=(0, 0, 0),
        #color=[1.0, 1.0, 0.9843, 1.0],  # Color temp. 5400 K
        intensity=1.,
        light_type='scene light')
    light2 = pv.Light(
        position=(-1000,-500,200),
        focal_point=(0, 0, 0),
        #color=[1.0, 1.0, 0.9843, 1.0],  # Color temp. 5400 K
        intensity=10,
        light_type='scene light')
    p.add_light(light1)
    #p.add_light(light2)



if 0:
    p.camera_position = 'xz'
    p.camera.azimuth = 160
    p.camera.elevation = 20




#p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,20))
toposurf = p.add_mesh(mapped_surf,texture=texture)


# big picture animation
p.camera_position =  [(-580.544, -4970.567, 3297.208), (-320.092, -61.711, -376.424), (0.005, 0.599, 0.801)]
p.camera_position =  [(-3218.235, -4414.535, 2936.848), (2.647, -244.782, -209.476), (0.226, 0.469, 0.854)]

# zoom animation
#p.camera_position =  [(-631.264, -407.928, 608.236), (-618.312, 380.175, 149.078), (-0.01, 0.504, 0.864)]tower

add_tower = True

if add_tower:
    tower = pv.read('/Users/rjl/B/pyvista/towers/OSVES3Dmodel.stl')
    xtower,ytower = -124.113, 46.9046
    scale = 0.3048  # convert from feet to meters

    itower = where(fgout.x <= xtower)[0].max()
    jtower = where(fgout.y <= ytower)[0].max()
    Btower = fgout.B[itower,jtower]
    print('Btower = %.1fm' % Btower)

    xtower_shift = (xtower - x_mean) * 111e3 * cos(y_mean*pi/180)
    ytower_shift = (ytower - y_mean) * 111e3
    ztower_shift = (Btower - 0.3) * warpfactor

    tower.points *= scale
    tower.points[:,0] += xtower_shift
    tower.points[:,1] += ytower_shift
    tower.points[:,2] *= warpfactor
    tower.points[:,2] += ztower_shift

    p.add_mesh(tower, color='r')

add_lighthouse = True

if add_lighthouse:
    lhouse = pv.read('/Users/rjl/B/pyvista/towers/lighthouse/scene.gltf')
    xlhouse,ylhouse = -124.1169, 46.8882
    scale = 0.3048*107 / 11.  # convert from feet to meters

    ilhouse = where(fgout.x <= xlhouse)[0].max()
    jlhouse = where(fgout.y <= ylhouse)[0].max()
    Blhouse = fgout.B[ilhouse,jlhouse]
    print('Blhouse = %.1fm' % Blhouse)

    xlhouse_shift = (xlhouse - x_mean) * 111e3 * cos(y_mean*pi/180)
    ylhouse_shift = (ylhouse - y_mean) * 111e3
    zlhouse_shift = (Blhouse - 0.3) * warpfactor

    def fix_lhouse(grid):

        # scale to desired size:
        grid.points *= scale

        # swap y and z to rotate to be vertical:
        y = grid.points[:,2].copy()
        z = grid.points[:,1].copy()
        grid.points[:,1] = y
        grid.points[:,2] = z

        # translate to desired location:
        grid.points[:,0] += xlhouse_shift
        grid.points[:,1] += ylhouse_shift
        grid.points[:,2] += zlhouse_shift

    def iter_grids(mb,fcn,level=0):
        if mb.keys() == [None]:
            grid = mb[0]
            print('Fixing grid at level %i with %i points' \
                    % (level,len(grid.points)))
            fcn(grid)
        else:
            for k in mb.keys():
                #print('Recursive call on level %i...' % (level+1))
                iter_grids(mb[k], fcn, level+1)

    iter_grids(lhouse, fix_lhouse)

    p.add_mesh(lhouse, color='r')

if 1:
    topoxyz.point_data['eta'] = eta.flatten(order='F')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    etamesh = p.add_mesh(etawarp,color='c')


def set_frameno(fgframeno):
    global etamesh, toposurf
    #fgframeno = 6*int(floor(tmin)) + 1
    fgframeno = int(round(fgframeno))
    fgout = fgout_grid.read_frame(fgframeno)
    thour, remainder = divmod(fgout.t, 3600)
    tmin, tsec = divmod(remainder, 60)
    tstr = '%s:%s:%s' % (str(int(thour)).zfill(2),
                      str(int(tmin)).zfill(2),
                      str(int(tsec)).zfill(2))
    print('Frame %i, t = %s after %s earthquake' % (fgframeno, tstr,event))


    if 0:
        # replot topo in case of co-seismic dz
        B = fgout.B
        B = flipZ(B)
        B = minimum(B, Bmax)
        topoxyz.point_data['B'] = B.flatten(order='F')
        topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)
        p.remove_actor(toposurf)
        toposurf = p.add_mesh(mapped_surf,texture=texture)


    #eta = where(fgout.h>0.1, fgout.eta, nan)
    eta = where(fgout.h>0.1, fgout.h + B0, nan)
    #eta = where(B0>0, fgout.h, eta)  # replace eta by h on original shore
    eta = flipZ(eta)

    if add_tower:
        print('At tower: h = %.2fm, eta = %.2fm, B = %.2f' \
                % (fgout.h[itower,jtower], fgout.eta[itower,jtower],
                   fgout.B[itower,jtower]))

    topoxyz.point_data['eta'] = eta.flatten(order='F')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)
    etamesh = p.add_mesh(etawarp,
                         scalars='eta',cmap=cmap_water,clim=(0,5))
                         #color='c')
                     #scalars='eta',cmap=geoplot.tsunami_colormap,clim=(-5,5))

    p.add_title('Westport at time %s after %s earthquake' % (tstr,event))

    if not make_animation:
        # round off entries in p.camera_position and print out, so user
        # can copy and paste into this script once good position is found:
        camera_position = list(p.camera_position)
        for i,a in enumerate(camera_position):
            b = []
            for j in range(len(a)):
                b.append(round(a[j],3))
            camera_position[i] = tuple(b)
        print('p.camera_position = ',camera_position)

    # to capture each frame after moving slider, uncomment this:
    #p.screenshot('PyVista_Frame%s.png' % str(fgframeno).zfill(4))



if not make_animation:
    #p.show_bounds(grid=True)
    p.add_slider_widget(set_frameno, [1,360], value=1, title='Frame',
                    pointa=(0.75,0.2), pointb=(0.95,0.2), color='k',
                    slider_width=0.02, tube_width=0.0025)

    p.show()

else:
    # animation:
    #framenos = range(1,362,3) # big picture
    framenos = range(1,182,3) # incomplete run
    #framenos = range(88,122,1) # for zoom

    for fgframeno in framenos:
        set_frameno(fgframeno)
        fname_png = '%s/PyVistaFrame%s.png' % (framedir, str(fgframeno).zfill(4))
        p.screenshot(fname_png)
        print('Created ',fname_png)

    p.close()

    anim = animation_tools.make_anim(framedir, fname_pattern='PyVistaFrame*.png')
    fname_mp4 = 'WestportN_animation_%s.mp4' % event
    animation_tools.make_mp4(anim, fname_mp4)
    print('Created ',fname_mp4)
