"""
Adapted from
    Hoquiam_2024/Runs/CSZ_L1_extended_pmel_SLR0111_BIG/pyvista_fgout2_animate.py

Need to include B0 for colormap
"""


from pylab import *
import os
import pyvista as pv
from clawpack.geoclaw import fgout_tools
from clawpack.visclaw import animation_tools
from clawpack.visclaw import colormaps, geoplot
from time import sleep

make_animation = True
fname_mp4 = 'Flattery_animation_CSZL1.mp4'   # used if make_animation

fgno = 3  # which fgout grid

outdir = '_output'
format = 'binary'  # format of fgout grid output

# where to save a png file for each frame, for constructing an animation:
framedir = '_frames'
os.system('mkdir -p %s' % framedir)
os.system('rm %s/*' % framedir)  # remove frames from previous version


#p = pv.Plotter(off_screen=make_animation, lighting='none')
#p = pv.Plotter(off_screen=make_animation, lighting='three lights')
p = pv.Plotter(off_screen=make_animation)

#p.window_size=(2300,1300)
p.window_size=(2800,1600)

if 0:
    light = pv.Light()
    light.set_direction_angle(30, 20)
    p.add_light(light)
    light = pv.Light()
    light.set_direction_angle(30, 0)
    p.add_light(light)

if 0:
    p.camera_position = 'xz'
    p.camera.azimuth = 180
    p.camera.elevation = 40


cpos =  [(3683.0523119722297, 6550.809109197466, 9805.097079591535),
 (5711.868885110077, -4237.337231554535, 1474.1260495398374),
 (0.024170198220371836, -0.6081729926498342, 0.7934364577767601)]

# Instantiate object for reading fgout frames:
fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format)
fgout_grid.read_fgout_grids_data(fgno)

# read solution eta and topo B at initial time (assumed pre-seismic)
fgout = fgout_grid.read_frame(1)

# convert x,y to meters, roughly:

x = (fgout.x - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
y = -(fgout.y - fgout.y[0]) * 111e3
z = array([0.])
X,Y,Z = meshgrid(x, y, z, indexing='ij')
topoxyz = pv.StructuredGrid(X,Y,Z)

warpfactor = 3  # amplification of elevations
Bmax = 500.  # truncate hills at this elevation

B0 = fgout.B  # save this as B0 for later plots of h+B0
B = flipud(B0)
B = minimum(B, Bmax)
topoxyz.point_data['B'] = B.flatten(order='F')
topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)


if 0:
    # plot topo alone:
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-100,100))
    p.view_xy()
    p.show(window_size=(1500,1500))

use_texture = False  # drape an image on the topo?  If False, color by elevation

if use_texture:
    texture = pv.read_texture('../Flattery_img.jpg')
    GE_extent = [-124.76, -124.55, 48.29, 48.41]
    #GE_extent = fgout.extent_edges  # must agree for texture mapping?

    x2,x1 = (asarray(GE_extent[:2]) - fgout.x[0]) * 111e3 * cos(fgout.y.mean()*pi/180)
    y2,y1 = (asarray(GE_extent[2:]) - fgout.y[0]) * 111e3

    o = (x1, y1, 0.)
    u = (x2, y1, 0.)
    v = (x1, y2, 0.)

    mapped_surf = topowarp.texture_map_to_plane(o, u, v)


if use_texture:
    toposurf = p.add_mesh(mapped_surf,texture=texture)
else:
    toposurf = p.add_mesh(topowarp,cmap='gist_earth',clim=(-0.8*Bmax, 0.8*Bmax))

if 0:
    #p.view_xy()
    cpos =  [(7525.908555910348, 27200.234921178228, 43445.067939385546),
     (6914.234289199344, -4076.481626872982, 1492.5276702953415),
     (-0.009641844655725214, -0.8016137124617913, 0.5977645780948045)]
    cpos = p.show(cpos=cpos,return_cpos=True)  # show topo alone for debugging
    print('cpos = ',cpos)

# initial sea surface eta:
# replace land surface by nan in eta so it only shows water:
eta = where(fgout.h>0.1, fgout.eta, nan)
eta = flipud(eta)

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
    print('Frame %i, t = %s after CSZ L1 event' % (fgframeno, tstr))
    if 1:
        # replot topo in case of co-seismic dz
        B = fgout.B
        B = flipud(B)
        B = minimum(B, Bmax)
        topoxyz.point_data['B'] = B.flatten(order='F')
        topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)
        p.remove_actor(toposurf)
        if use_texture:
            toposurf = p.add_mesh(mapped_surf,texture=texture)
        else:
            toposurf = p.add_mesh(topowarp,cmap='gist_earth',
                                  clim=(-0.8*Bmax, 0.8*Bmax))

    # plot eta relative to deformed topo:
    eta = where(fgout.h>0.1, fgout.eta, nan)

    # plot eta relative to original topo B0:
    # (so colors still show eta relative to original shore elevation)
    #eta = where(fgout.h>0.1, fgout.h + B0, nan)

    #eta = where(fgout.h>0.1, fgout.eta - fgout.B + B0, nan)

    #eta = where(B0>0, fgout.h, eta)  # replace eta by h on original shore

    eta = flipud(eta)
    topoxyz.point_data['eta'] = eta.flatten(order='F')
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
    p.remove_actor(etamesh)

    etamesh = p.add_mesh(etawarp,scalars='eta',
                         cmap=geoplot.tsunami_colormap,clim=(-5,5))

    p.add_title('Time %s after CSZ L1 event' % tstr)

    if not make_animation and False:
        # NOW RETURNED BY show()
        # print camera position so that this can be copied and pasted
        # into this script after adjusting (and then sliding frameno)
        print('p.camera_position = ', p.camera_position)

    # to capture each frame after moving slider, uncomment this:
    #p.screenshot('PyVista_Frame%s.png' % str(fgframeno).zfill(4))


if not make_animation:
    p.add_slider_widget(set_frameno, [1,182], value=1, title='Frame',
                    pointa=(0.75,0.85), pointb=(0.95,0.85), color='gray',
                    slider_width=0.02, tube_width=0.005)

    cpos = p.show(cpos=cpos,return_cpos=True)  # show topo alone for debugging
    print('Final camera_position: \ncpos = ',cpos)

else:

    p.camera_position = cpos
    for fgframeno in range(1,182,1):
        set_frameno(fgframeno)
        fname_png = '%s/PyVistaFrame%s.png' % (framedir, str(fgframeno).zfill(4))
        p.screenshot(fname_png)
        print('Created ',fname_png)

    p.close()

    anim = animation_tools.make_anim(framedir, fname_pattern='PyVistaFrame*.png')
    animation_tools.make_mp4(anim, fname_mp4)
    print('Created ',fname_mp4)
