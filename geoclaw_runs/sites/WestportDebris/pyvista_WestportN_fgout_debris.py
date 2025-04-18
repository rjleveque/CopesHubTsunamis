"""
Use PyVista to plot the tsunami in Westport with debris tracking
"""

from pylab import *
import os
import pyvista as pv
from clawpack.geoclaw import fgout_tools
import debris_tools
from clawpack.visclaw import animation_tools
from WestportN_debris import compute_debris_paths

global etamesh, debris_spheres, debris_path_t

rundir = os.getcwd()
event = os.path.split(rundir)[-1]

fgno = 1  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output
qmap = 'geoclaw'  # defines mapping for fgout variables


# where to save a png file for each frame, for constructing an animation:
framedir = '_frames'
os.system('mkdir -p %s' % framedir)
os.system('rm %s/*' % framedir)  # remove frames from previous version

#tfinal = 4780.  # create animation up to this time (or end of computation)
tfinal = 1.5*3600.

# compute debris paths based on fgout, based on buildings in this extent:
buildings_extent = [-124.129, -124.112, 46.8805, 46.8905]
#buildings_extent = None # to use all buildings

debris_paths, dbnos = compute_debris_paths(tfinal, buildings_extent)

#debris_spheres = len(debris_paths)*[' ']  # for plot actors
debris_spheres = {}  # for plot actors

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format, qmap)
#fgout_grid.read_fgout_grids_data_pre511(fgno)
try:
    fgout_grid.read_fgout_grids_data(fgno)
except:
    fgout_grid.read_fgout_grids_data_pre511(fgno)
    fgout_grid.times = linspace(fgout_grid.tstart, fgout_grid.tend,
                                fgout_grid.nout)
fgout = fgout_grid.read_frame(1)

# convert x,y to meters, roughly:

x = (fgout.x - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
y = -(fgout.y - fgout.y[0]) * 111e3
z = array([0.])
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

B = fgout.B
B = flipud(B)
Bmax = 25.
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='F')

warpfactor = 1  # amplification of elevations

topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

GE_file = '../GE_WestportN.jpg'
GE_extent = fgout.extent_edges
texture = pv.read_texture(GE_file)
x2,x1 = (asarray(GE_extent[:2]) - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
y2,y1 = (asarray(GE_extent[2:]) - fgout.y[0]) * 111e3

#x1,x2,y1,y2 = GE_extent

o = (x1, y1, 0.)
u = (x2, y1, 0.)
v = (x1, y2, 0.)

mapped_surf = topowarp.texture_map_to_plane(o, u, v)

if 0:
    # plot topo flat:
    p = pv.Plotter()
    p.add_mesh(topoxyz,cmap='gist_earth',clim=(-20,10))
    p.view_xy()
    p.show(window_size=(1500,1500))

if 0:
    # plot topo alone:
    p = pv.Plotter()
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,20))
    p.camera_position = 'xz'
    p.camera.azimuth = 180
    p.camera.elevation = 20
    p.show(window_size=(1500,1500))

# replace land surface by nan in eta so it only shows water:
eta = where(fgout.h>0.1, fgout.eta, nan)
eta = flipud(eta)

make_html = False  # True to make html version of topo (requires trame)

make_animation = True

# select frames for fgout grid 5 every 20 seconds up to tfinal:
fgframes = [n+1 for n in range(len(fgout_grid.times)) \
               if ((fgout_grid.times[n] <= tfinal) and \
                  (mod(fgout_grid.times[n], 20) == 0))]


# =================

p = pv.Plotter(off_screen=make_animation, lighting='three lights')

#p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,20))
p.add_mesh(mapped_surf,texture=texture)

topoxyz.point_data['eta'] = eta.flatten(order='F')
etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
etamesh = p.add_mesh(etawarp,color='c')

if 0:
    p.camera_position = 'xz'
    p.camera.azimuth = 120
    p.camera.elevation = 40


#p.camera_position =  [(-3218.235, -4414.535, 2936.848), (2.647, -244.782, -209.476), (0.226, 0.469, 0.854)]

p.camera_position =  [(3833.0778858778876, -62.37740956595462, 1546.5793815421837),
 (2709.0180605962937, -2255.950695273429, 9.55961395021373),
 (-0.21428121046481832, -0.4844864786159533, 0.8481488164703442)]

p.camera_position =  [(4021.4254960355893, 3013.3560518492995, 3063.9763892981005),
 (2725.9018709724173, -2103.254422467251, -239.49720889260217),
 (-0.06807467510541654, -0.5289284527556626, 0.845931752846997)]

p.window_size=(2300,1300)

debris_path_t = {}
for dbno in dbnos:
    debris_path_t[dbno] = None  # to indicate this debris not yet drawn

def set_frameno(fgframeno):
    global etamesh, debris_spheres, debris_path_t
    #fgframeno = 6*int(floor(tmin)) + 1
    fgframeno = int(round(fgframeno))
    fgout = fgout_grid.read_frame(fgframeno)
    thour, remainder = divmod(fgout.t, 3600)
    tmin, tsec = divmod(remainder, 60)
    tstr = '%s:%s:%s' % (str(int(thour)).zfill(2),
                      str(int(tmin)).zfill(2),
                      str(int(tsec)).zfill(2))
    print('Frame %i, t = %s' % (fgframeno, tstr))

    eta = where(fgout.h>0.1, fgout.eta, nan)
    eta = flipud(eta)
    topoxyz.point_data['eta'] = eta.flatten(order='F')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)
    etamesh = p.add_mesh(etawarp,color='c')

    p.add_title('Time %s after %s earthquake' % (tstr,event))

    # add debris particles at this time:
    #xdebris,ydebris = debris_tools.get_debris_xy(fgout.t, debris_paths, dbnos)

    new_debris_path_t = debris_tools.get_debris_path_t(fgout.t, debris_paths, dbnos)
    # returns dict with key dbno and value (tdebris,xdebris,ydebris,udebris,vdebris)
    # with tdebris == t (or nearest time the path exists)


    # make function that returns eta(x,y) at this time:
    eta_fcn = fgout_tools.make_fgout_fcn_xy(fgout, 'eta')
    B_fcn = fgout_tools.make_fgout_fcn_xy(fgout, 'B')


    #for k in range(len(xdebris)):
    #for k in range(0, len(xdebris), 1):
    for dbno in dbnos:
        #import pdb; pdb.set_trace()
        if (debris_path_t[dbno] is None):
            redraw_debris = True
        else:
            redraw_debris = (abs(new_debris_path_t[dbno][1:5] \
                                   - debris_path_t[dbno][1:5]).max() > 1e-5)
        if redraw_debris:
            # debris needs to be drawn, otherwise leave this artist alone
            xdebris,ydebris = new_debris_path_t[dbno][1:3]
            xd = -(xdebris- fgout.x[-1]) * 111e3 * cos(fgout.y.mean()*pi/180)
            yd = -(ydebris- fgout.y[0]) * 111e3
            boxl = sphere_radius = 10
            zd = warpfactor * max(eta_fcn(xdebris,ydebris),
                     B_fcn(xdebris,ydebris)+sphere_radius/warpfactor)
            #debris_sphere = pv.Sphere(radius=sphere_radius, center=(xd,yd,zd))
            debris_sphere = pv.Box((xd-boxl,xd+boxl,yd-boxl,yd+boxl,zd-boxl,zd+boxl))
            if debris_path_t[dbno] is not None:
                # remove debris from old location:
                p.remove_actor(debris_spheres[dbno])
            # draw in new location:
            debris_spheres[dbno] = p.add_mesh(debris_sphere, color='r')
            #print('+++ dbno = %i, redrawn at x = %.5f, y = %.5f' % (dbno,xd,yd))
        else:
            #print('+++ dbno = %i not redrawn' % dbno)
            pass
    debris_path_t = new_debris_path_t.copy()


    if not make_animation:
        # print camera position so that this can be copied and pasted
        # into this script after adjusting (and then sliding frameno)
        print('p.camera_position = ', p.camera_position)


    # to capture each frame after moving slider, uncomment this:
    #p.screenshot('PyVista_Nuu_Frame%s.png' % str(fgframeno).zfill(4))

#p.add_slider_widget(set_frameno, [0,20], value=0, title='Time (minutes)')
#p.add_slider_widget(set_frameno, [0,4], value=0, title='Time (minutes)')


if make_html:
    frameno_html = 62
    #frameno_html = 1
    set_frameno(frameno_html)
    fname_html = 'Nuu_frame%s.html' % str(frameno_html).zfill(4)
    p.export_html(fname_html)
    print('Created ',fname_html)


elif not make_animation:
    frames_range = [fgframes[0], fgframes[-1]]
    p.add_slider_widget(set_frameno, frames_range, value=1, title='Frameno',
                    pointa=(0.4,0.85), pointb=(0.9,0.85), color='blue',
                    slider_width=0.02, tube_width=0.005)

    p.show()

else:

    #for fgframeno in range(1,121,1):
    #for fgframeno in range(1,74,1):
    #for fgframeno in fgframes:
    for fgframeno in range(1,300,5):
        set_frameno(fgframeno)
        fname_png = '%s/PyVistaFrame%s.png' % (framedir, str(fgframeno).zfill(4))
        p.screenshot(fname_png)
        print('Created ',fname_png)

    p.close()

    anim = animation_tools.make_anim(framedir, fname_pattern='PyVistaFrame*.png')
    fname_mp4 = 'WestportN_debris_animation_%s.mp4' % event
    animation_tools.make_mp4(anim, fname_mp4)
    print('Created ',fname_mp4)
