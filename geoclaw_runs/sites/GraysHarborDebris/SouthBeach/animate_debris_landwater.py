
from pylab import *
import matplotlib.animation as animation
from matplotlib import colors
from scipy.interpolate import RegularGridInterpolator

from clawpack.visclaw import animation_tools

#from clawpack.geoclaw import fgout_tools
from geoclaw_debris import fgout_tools # modified version

from track_debris_points_landwater import compute_debris_paths


try:
    from geoclaw_debris import debris_tools
except:
    msg = '*** You need to add the directory geoclaw_debris/src to your ' \
          'Python path, from this repo:\n' \
          '     https://github.com/rjleveque/geoclaw_debris'
    raise ImportError(msg)
        

if 0:
    graphics_dir = '.'
    fgout_extent = [-124.165, -123.97, 46.73, 46.92] # SouthBeach
    bgimage = [imread(f'{graphics_dir}/SouthBeach_img.jpg'), fgout_extent]
else:
    bgimage = None

fgno = 1

event = 'BL10S'
#event = 'BR10M'
outdir = f'_output_{event}'

output_format = 'binary'  # format of fgout grid output


# load all fgout frames:
# Note this requires a modified fgout_tools.py since this function is not
# currently in the GeoClaw v5.14.0 version

qois = ['h', 'u', 'v', 'B']  # arrays to return
x,y,t,qois,qoi_arrays = fgout_tools.load_fgout_qois(fgno, outdir,
                                                     output_format, qois)

X,Y = meshgrid(x,y,indexing='ij')  # should agree with fgout X and Y

if 0:
    # List of times to use for making debris paths (and animation):
    #db_i = range(20,60,10)
    db_i = range(0,180,2)
    db_times = t[db_i]
else:
    # use all frames:
    db_times = t

print(f'Will move debris using these discrete times: \n{db_times}')


# make interpolating functions, e.g. qoi_fcns['h'] for depth at (x,y,t)
qoi_fcns = {}
for qoi in qois:
    qoi_fcns[qoi] =fgout_tools.make_qoi_fcn_xyt(x,y,t,qoi_arrays[qoi])


fgout_extent = [x.min(), x.max(), y.min(), y.max()]

debris_paths, dbnos_water, dbnos_land = compute_debris_paths(db_times,
                                                             qoi_fcns)

if 1:
    # make plots and animation:

    qoi = 'h'  # 'h' or 'speed' (currently set up only for 'h')

    B0 = qoi_fcns['B'](X, Y, db_times[0])
    h0 = qoi_fcns['h'](X, Y, db_times[0])
    t0 = db_times[0]

    Bshore = 0.  # FIX

    ylat = Y.mean()  # for aspect ratio of plots

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

        im = imshow(flipud(h0.T), extent=fgout_extent,
                    cmap=cmap_depth, norm=norm_depth)
        cb = colorbar(im, extend='max', shrink=0.7)
        cb.set_label('meters')
        contour(X, Y, B0, [Bshore], colors='g', linewidths=0.5)

        gca().set_aspect(1./cos(ylat*pi/180.))
        ticklabel_format(useOffset=False)
        xticks(rotation=20)

        title_text= title(f'Water depth {t0/60:.2f} minutes after quake')


        points_water, = debris_tools.plot_debris(t0, debris_paths,
                                                dbnos_water, ax,
                                                color='b',
                                                marker='s',markersize=0.5)
        
        points_land, = debris_tools.plot_debris(t0, debris_paths,
                                                dbnos_land, ax,
                                                color='r',
                                                marker='s',markersize=0.5)



    def update(t):
        """
        Update the data in the plot for a new time
        """

        h = qoi_fcns['h'](X, Y, t)
        im.set_data(flipud(h.T))

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
    #anim_times = db_times  # same set of times used for tracking
    anim_times = db_times[::2]
    anim = animation.FuncAnimation(fig, update,
                          frames=anim_times, interval=200, blit=False)

    fname_mp4 = f'landwater_tracers_{event}_new_xyt.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
