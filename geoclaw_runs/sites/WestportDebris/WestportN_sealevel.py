from pylab import *
import pyvista as pv
from clawpack.geoclaw import topotools

topo = topotools.Topography('/Users/rjl/topo/topofiles/GH_tiles_2021_13s.asc')

extent = [-124.14, -124.07, 46.88, 46.92]  # WestportN

topo = topo.crop(extent)

#topo.plot(limits=(-25,25))

z = array([0.])
x_mean = topo.x.mean()
y_mean = topo.y.mean()
x = (topo.x - x_mean) * 111e3 * cos(y_mean*pi/180)
y = -(topo.y - y_mean) * 111e3
#x = (topo.x - topo.x[0]) * 111e3 * cos(y_mean*pi/180)
#y = -(topo.y - topo.y[0]) * 111e3
print('topo coordinates in meters: ')
print('     xmin = %.1fm, xmax = %.1fm' % (x.min(),x.max()))
print('     ymin = %.1fm, ymax = %.1fm' % (y.min(),y.max()))
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

#B = topo.Z
B = flipud(topo.Z)
#B = fliplr(B)
Bmax = 50.
#B = where(B<Bmax,B,nan)
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='C')

warpfactor = 1 #3 #15  # amplification of elevations

topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

if 0:
    # plot topo alone:
    p = pv.Plotter()
    #p.add_mesh(topoxyz,cmap='gist_earth',clim=(-5,5))
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,10))
    #p.add_mesh(topowarp,cmap='terrain',clim=(-20,10))
    p.camera_position =  [(1305.71, -7697.221, 4709.439), (713.629, -2677.6, 910.373), (-0.001, 0.603, 0.797)]
    p.show(window_size=(1500,1500))



global etamesh

make_snapshots = False
#p = pv.Plotter()
p = pv.Plotter(off_screen=make_snapshots)
#p.add_mesh(topowarp,cmap='gist_earth',clim=(-5,5))

# Add GE image as texture:

GE_file = 'GE_WestportN.jpg'
GE_extent = extent

texture = pv.read_texture(GE_file)

x1,x2 = (asarray(GE_extent[:2]) - x_mean) * 111e3 * cos(topo.y.mean()*pi/180)
y1,y2 = (asarray(GE_extent[2:]) - y_mean) * 111e3

o = (x1, y1, 0.)
u = (x2, y1, 0.)
v = (x1, y2, 0.)

mapped_surf = topowarp.texture_map_to_plane(o, u, v)
p.add_mesh(mapped_surf,texture=texture)


add_tower = True

if add_tower and 0:
    # use OSVES tower:
    tower = pv.read('/Users/rjl/B/pyvista/towers/OSVES3Dmodel.stl')
    xtower,ytower = -124.113, 46.9046
    scale = 0.3048  # convert from feet to meters

    itower = where(topo.x <= xtower)[0].max()
    jtower = where(topo.y <= ytower)[0].max()
    Btower = topo.Z[jtower,itower]
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

else:
    # generic tower
    tower = pv.read('/Users/rjl/B/pyvista/towers/tower_house_design/scene.gltf')
    scale = 40

    xtower,ytower = -124.113, 46.9046

    i = where(topo.x <= xtower)[0].max()
    j = where(topo.y <= ytower)[0].max()
    Btower = topo.Z[j,i]
    print('Btower = %.1fm' % Btower)

    xtower_shift = (xtower - x_mean) * 111e3 * cos(y_mean*pi/180)
    ytower_shift = (ytower - y_mean) * 111e3
    ztower_shift = (Btower - 0.3) * warpfactor

    #import pdb; pdb.set_trace()

    global xmin,ymin,zmin, xmax,ymax,zmax


    def find_limits(grid):
        global xmin,ymin,zmin, xmax,ymax,zmax
        xmin = min(xmin, grid.points[:,0].min())
        xmax = max(xmax, grid.points[:,0].max())
        ymin = min(ymin, grid.points[:,1].min())
        ymax = max(ymax, grid.points[:,1].max())
        zmin = min(zmin, grid.points[:,2].min())
        zmax = max(zmax, grid.points[:,2].max())

    def fix_tower(grid):

        # scale to desired size:
        grid.points *= scale

        # swap y and z to rotate to be vertical:
        y = grid.points[:,2].copy()
        z = grid.points[:,1].copy()
        grid.points[:,1] = y
        grid.points[:,2] = z

        # translate to desired location:
        grid.points[:,0] += xtower_shift
        grid.points[:,1] += ytower_shift
        grid.points[:,2] += ztower_shift

    def iter_grids(mb,fcn,level=0):
        if mb.keys() == [None]:
            grid = mb[0]
            if 0:
                print('Fixing grid at level %i with %i points' \
                        % (level,len(grid.points)))
            fcn(grid)
        else:
            for k in mb.keys():
                #print('Recursive call on level %i...' % (level+1))
                iter_grids(mb[k], fcn, level+1)

    xmin = inf; ymin = inf; zmin = inf
    xmax = -inf; ymax = -inf; zmax = -inf
    iter_grids(tower, find_limits)
    print('tower xmin = %10.3f,    xmax = %10.3f,    x-length = %10.3f' \
            % (xmin,xmax, xmax-xmin))
    print('tower ymin = %10.3f,    ymax = %10.3f,    y-length = %10.3f' \
            % (zmin,zmax, ymax-ymin))
    print('tower zmin = %10.3f,    zmax = %10.3f,    z-height = %10.3f' \
            % (ymin,ymax, zmax-zmin))

    iter_grids(tower, fix_tower)

    xmin = inf; ymin = inf; zmin = inf
    xmax = -inf; ymax = -inf; zmax = -inf
    iter_grids(tower, find_limits)
    print('after scaling, shifting, converting to meters:')
    print('tower xmin = %10.3f,    xmax = %10.3f,    x-length = %10.3f' \
            % (xmin,xmax, xmax-xmin))
    print('tower ymin = %10.3f,    ymax = %10.3f,    y-length = %10.3f' \
            % (ymin,ymax, ymax-ymin))
    print('tower zmin = %10.3f,    zmax = %10.3f,    z-height = %10.3f' \
            % (zmin,zmax, zmax-zmin))
    xtower2 = x_mean + 0.5*(xmin+xmax) / (111e3*cos(y_mean*pi/180))
    ytower2 = y_mean + 0.5*(ymin+ymax) / 111e3
    print('xtower  = %11.5f,  ytower  = %11.5f' % (xtower,ytower))
    print('xtower2 = %11.5f,  ytower2 = %11.5f' % (xtower2,ytower2))

    p.add_mesh(tower, color='y')



add_lighthouse = True

if add_lighthouse:
    lhouse = pv.read('/Users/rjl/B/pyvista/towers/lighthouse/scene.gltf')
    xlhouse,ylhouse = -124.1169, 46.8882
    scale = 0.3048*107 / 11.  # convert from feet to meters

    ilhouse = where(topo.x <= xlhouse)[0].max()
    jlhouse = where(topo.y <= ylhouse)[0].max()
    Blhouse = topo.Z[jlhouse,ilhouse]
    print('Blhouse = %.1fm' % Blhouse)

    xlhouse_shift = (xlhouse - x_mean) * 111e3 * cos(y_mean*pi/180)
    ylhouse_shift = (ylhouse - y_mean) * 111e3
    zlhouse_shift = (Blhouse - 0.3) * warpfactor

    def find_limits(grid):
        global xmin,ymin,zmin, xmax,ymax,zmax
        xmin = min(xmin, grid.points[:,0].min())
        xmax = max(xmax, grid.points[:,0].max())
        ymin = min(ymin, grid.points[:,1].min())
        ymax = max(ymax, grid.points[:,1].max())
        zmin = min(zmin, grid.points[:,2].min())
        zmax = max(zmax, grid.points[:,2].max())

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


    if 0:
        xmin = inf; ymin = inf; zmin = inf
        xmax = -inf; ymax = -inf; zmax = -inf
        iter_grids(tower, find_limits)
        print('tower xmin = %10.3f,    xmax = %10.3f,    x-length = %10.3f' \
                % (xmin,xmax, xmax-xmin))
        print('tower ymin = %10.3f,    ymax = %10.3f,    y-length = %10.3f' \
                % (zmin,zmax, ymax-ymin))
        print('tower zmin = %10.3f,    zmax = %10.3f,    z-height = %10.3f' \
                % (ymin,ymax, zmax-zmin))

        iter_grids(tower, fix_tower)

        xmin = inf; ymin = inf; zmin = inf
        xmax = -inf; ymax = -inf; zmax = -inf
        iter_grids(tower, find_limits)
        print('after scaling, shifting, converting to meters:')
        print('tower xmin = %10.3f,    xmax = %10.3f,    x-length = %10.3f' \
                % (xmin,xmax, xmax-xmin))
        print('tower ymin = %10.3f,    ymax = %10.3f,    y-length = %10.3f' \
                % (ymin,ymax, ymax-ymin))
        print('tower zmin = %10.3f,    zmax = %10.3f,    z-height = %10.3f' \
                % (zmin,zmax, zmax-zmin))
        xtower2 = x_mean + 0.5*(xmin+xmax) / (111e3*cos(y_mean*pi/180))
        ytower2 = y_mean + 0.5*(ymin+ymax) / 111e3
        print('xtower  = %11.5f,  ytower  = %11.5f' % (xtower,ytower))
        print('xtower2 = %11.5f,  ytower2 = %11.5f' % (xtower2,ytower2))



sea_level = 0.
eta = where(B < sea_level, sea_level, nan)
topoxyz.point_data['eta'] = eta.flatten(order='C')
etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
etamesh = p.add_mesh(etawarp,color='c')


p.camera_position =  [(-1359.499, -3780.495, 1047.518), (-524.271, 50.108, -434.935), (0.042, 0.353, 0.935)]

p.camera_position =  [(-3218.235, -4414.535, 2936.848), (2.647, -244.782, -209.476), (0.226, 0.469, 0.854)]

def set_sea_level(sea_level):
    global etamesh
    #p.clear() # causes seg fault
    eta = where(B < sea_level, sea_level, nan)
    topoxyz.point_data['eta'] = eta.flatten(order='C')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)
    etamesh = p.add_mesh(etawarp,color='c')
    p.add_title('Westport at MHW after %.2f m subsidence (or sea level rise)' \
                % sea_level)

    if 0:
        # test Rectangle
        xd =  (-124.12 - x_mean) * 111e3 * cos(y_mean*pi/180)
        yd =  (46.89 - y_mean) * 111e3
        zd = 5
        boxl = 100
        pA = [xd-boxl, yd-boxl, zd+warpfactor]
        pB = [xd-boxl, yd+boxl, zd+warpfactor]
        pC = [xd+boxl, yd+boxl, zd+warpfactor]
        debris_square = pv.Rectangle([pA,pB,pC])
        p.add_mesh(debris_square, color='r')

    if 1:
        # round off entries in p.camera_position and print out, so user
        # can copy and paste into this script once good position is found:
        camera_position = list(p.camera_position)
        for i,a in enumerate(camera_position):
            b = []
            for j in range(len(a)):
                b.append(round(a[j],3))
            camera_position[i] = tuple(b)
        print('p.camera_position = ',camera_position)

if not make_snapshots:
    p.add_slider_widget(set_sea_level, [-5,20], value=0, title='Sea Level',
                        pointa=(0.75,0.1), pointb=(0.95,0.1), color='k',
                        slider_width=0.02, tube_width=0.0025)

    light = pv.Light(position=(20,20,50))
    light.positional = True
    #p.add_light(light)
    p.show(window_size=(2500,1500))
    #p.export_html('OSVES_interactive.html')
else:
    p.window_size = (2500,1500)
    for slr in [0,1,2,3]:
        set_sea_level(slr)
        p.add_title('Westport at MHW after %.2f m subsidence (or sea level rise)' % slr)
        fname_png = 'Hoquiam_mhw%im.png' % slr
        p.screenshot(fname_png)
        print('Created ',fname_png)
    p.close()
