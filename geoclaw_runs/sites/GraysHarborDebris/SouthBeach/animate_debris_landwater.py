
from pylab import *
import matplotlib.animation as animation
from matplotlib import colors

from clawpack.visclaw import animation_tools
from clawpack.geoclaw import fgout_tools
from track_debris_points_landwater import compute_debris_paths
from geoclaw_debris import debris_tools


fgno = 1

if 0:
    graphics_dir = '/Users/rjl/git/WestportMaritime/graphics/'
    fgout_extent = [-124.16, -124.08, 46.885, 46.92]  # Westport
    bgimage = [imread(graphics_dir+'fgout11CT.png'), fgout_extent]
else:
    bgimage = None

qoi = 'h'  # 'h' or 'speed'

#event = 'BL10S'
event = 'BR10M'
outdir = f'_output_{event}'

output_format = 'binary'  # format of fgout grid output
qmap = 'geoclaw'  # defines mapping for fgout variables


# List of frames to use for making debris paths and animation:
fgframes = range(20,40,10)
#fgframes = range(1,180,2)

fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, output_format, qmap)
fgout_grid.read_fgout_grids_data(fgno)

# read all specified fgout frames (could modify to read single netCDF version):
# TODO

fgout0 = fgout_grid.read_frame(fgframes[0])
t0 = fgout0.t
fgout_extent = fgout0.extent_edges

debris_paths, dbnos_water, dbnos_land = compute_debris_paths(fgframes,
                                                fgno, outdir,
                                                output_format,
                                                extent=None)

if 1:
    # make plots and animation:
    fgout = fgout0
    B0 = fgout0.B  # assumes no subsidence
    Bshore = 0.  # FIX

    ylat = fgout.Y.mean()  # for aspect ratio of plots

    fig,ax = subplots(figsize=(7,7))
    if bgimage:
        ax.imshow(bgimage[0],extent=bgimage[1])

    if qoi=='h':
        ## Plot maximum flow depth

        bounds_depth = array([1e-6,0.5,1.0,1.5,2,4.0,6.0])
        bounds_depth = array([1e-6,1,2,3])

        a = 0.4 # alpha (transparency,  1=opaque)
        cmap_depth = colors.ListedColormap([[.7,.7,1],[.5,.5,1],[0,0,1],
                         [1,.7,.7], [1,.4,.4], [1,0,0]])
        cmap_depth = colors.ListedColormap([[.5,.5,1,a],[.3,.3,1,a],[0,0,1,a]])
        cmap_depth = colors.ListedColormap([[.3,.3,1,a],[0,0,1,a],[.5,0,.5,a]])

        # Set color for value exceeding top of range to purple:
        #cmap_depth.set_over(color=[1,0,1])


        if bgimage:
            # Set color to transparent where no water:
            cmap_depth.set_under(color=[1,1,1,0])
        else:
            # Set color to light green where no water:
            cmap_depth.set_under(color=[.7,1,.7,a])

        norm_depth = colors.BoundaryNorm(bounds_depth, cmap_depth.N)

        im = imshow(flipud(fgout0.h.T), extent=fgout_extent,
                    cmap=cmap_depth, norm=norm_depth)
        cb = colorbar(im, extend='max', shrink=0.7)
        cb.set_label('meters')
        contour(fgout.X, fgout.Y, B0, [Bshore], colors='g', linewidths=0.5)

        gca().set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20)

        t = fgout.t
        title_text= title(f'Water depth {t/60:.2f} minutes after quake')


        points_water, = debris_tools.plot_debris(fgout.t, debris_paths,
                                                dbnos_water, ax,
                                                color='b',
                                                marker='s',markersize=0.5)
        
        points_land, = debris_tools.plot_debris(fgout.t, debris_paths,
                                                dbnos_land, ax,
                                                color='r',
                                                marker='s',markersize=0.5)



    def update(fgframe):
        fgout = fgout_grid.read_frame(fgframe)
        t = fgout.t
        im.set_data(flipud(fgout.h.T))

        xd = []
        yd = []
        for dbno in dbnos_water:
            db = debris_paths[dbno]
            try:
                j = where(abs(db[:,0]-t) < 1e-6)[0].max()
            except:
                print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
                j = -1
            if j > -1:
                xd.append(db[j,1])
                yd.append(db[j,2])
        points_water.set_data(xd, yd)

        xd = []
        yd = []
        for dbno in dbnos_land:
            db = debris_paths[dbno]
            try:
                j = where(abs(db[:,0]-t) < 1e-6)[0].max()
            except:
                print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
                j = -1
            if j > -1:
                xd.append(db[j,1])
                yd.append(db[j,2])
        points_land.set_data(xd, yd)

        title_text.set_text(f'Water depth {(t/60):.2f} minutes after quake')

    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                          frames=fgframes, interval=200, blit=False)

    fname_mp4 = f'landwater_tracers_{event}_new.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
