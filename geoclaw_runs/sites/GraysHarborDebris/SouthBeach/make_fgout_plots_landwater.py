
from pylab import *
import postproc_fgout as P
from clawpack.visclaw import animation_tools
import matplotlib.animation as animation
from matplotlib import colors



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

format = 'binary'  # format of fgout grid output

# Create ClawPlotData object used for reading in fgout frames:
plotdata = P.make_plotdata(fgno, outdir, format)

# List of frames to use for making debris paths and animation:
#fgframes = range(20,40,5)
fgframes = range(1,180,2)

# Deterime time t0 of first fgout frame, to initialize particles
frameno0 = fgframes[0]
fgout0 = P.read_fgout_frame(fgno,frameno0,plotdata)
t0 = fgout0.t


B0 = fgout0.B
if t0 == 0.:
    Bshore = 0.
else:
    Bshore = -2.  # approximate subsidence for land vs water

Bshore = 0.  # NEED TO FIX SUBSIDENCE

fgout_extent = fgout0.extent

# make a B(x,y) function that can be used to determine land vs. water:
B0_fcn_xy = P.make_fgout_fcn_xy(fgout0, 'B')

# Initialize debris_paths dictionary and set
# initial debris particle locations (x0,y0) at time t0.
# Require a list of dbnos and each array
#     debris_paths[dbno]
# in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

debris_paths = {}
dbnos_water = []
dbnos_land = []
grounding_depth = {}
drag_factor = {}

# set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
u0 = 0.
v0 = 0.
if 0:
    x1 = fgout_extent[0] + 0.001
    x2 = fgout_extent[1] - 0.0005
    y1 = fgout_extent[2] + 0.001
    y2 = fgout_extent[3] - 0.0005
x1,x2,y1,y2 = [-124.14, -124.07, 46.74, 46.92]

dxd = 0.002
dyd = dxd * cos(47*pi/180)
xd = arange(x1, x2, dxd)
yd = arange(y1, y2, dyd)
mxd = len(xd)
myd = len(yd)

for i in range(len(xd)):
    x0 = xd[i]
    for j in range(len(yd)):
        y0 = yd[j]
        dbno = myd*i + j
        db = array([[t0, x0, y0, u0, v0]])
        debris_paths[dbno] = db

        if B0_fcn_xy(x0,y0) < Bshore:
            dbnos_water.append(dbno)
            grounding_depth[dbno] = 0.
            drag_factor[dbno] = None
        else:
            dbnos_land.append(dbno)
            grounding_depth[dbno] = 0.5
            drag_factor[dbno] = None #0.1

dbnos = dbnos_water + dbnos_land
print(f'Created {len(dbnos)} initial debris particles')
print(f'      {len(dbnos_water)} in the water and {len(dbnos_land)} onshore')


# Compute debris path for each particle by using all the fgout frames
# in the list fgframes (first frame should be frameno0 used to set t0 above):

debris_paths = P.make_debris_paths(fgno, fgframes, plotdata,
                                   debris_paths, dbnos, drag_factor,
                                   grounding_depth)


if 0:
    def plot_frame(fgout):
        fig,ax = P.plot_fgout_frame(fgout, qoi, bgimage)
        P.plot_debris(fgout.t, debris_paths, dbnos_water, ax, color='b',
                    marker='.',markersize=2)
        P.plot_debris(fgout.t, debris_paths, dbnos_land, ax, color='r',
                    marker='.',markersize=2)
        ax.set_xlim(fgout_extent[:2])
        ax.set_ylim(fgout_extent[2:])
        ax.contour(fgout.X, fgout.Y, B0, [Bshore], colors='g', linewidths=0.5)
        ax.plot([-124.11537],[46.907685],color='g',marker='*',markersize=8)
        return fig,ax

    # Plot initial frame and debris locations:
    fig,ax = plot_frame(fgout0)

if 0:
    # Plot full debris paths on the previous plot:
    for k in range(len(dbnos_water)):
        debris_path = debris_paths[dbnos_water[k]]
        ax.plot(debris_path[:,1], debris_path[:,2], color='k', linewidth=2)
    for k in range(len(dbnos_land)):
        debris_path = debris_paths[dbnos_land[k]]
        ax.plot(debris_path[:,1], debris_path[:,2], color='r', linewidth=1)

    if 0:
        fname = 'debris_paths.png'
        savefig(fname)
        print('Created ',fname)

if 0:
    # make mp4 animation:
    #fname_mp4 = '%s_%s.mp4' % (plotdata.file_prefix,qoi)
    fname_mp4 = 'test.mp4'
    fps = 5 # frames per second
    P.make_mp4(fgno, fgframes, plotdata, plot_frame, fname_mp4, fps)

if 0:
    for frameno in fgframes:
        fgout = P.read_fgout_frame(fgno,frameno,plotdata)
        t = fgout.t
        fig,ax = plot_frame(fgout)
        fname = '_plots/fgout%s.png' % str(frameno).zfill(4)
        fig.savefig(fname, dpi=60)
        print('Created ',fname)
        close(fig)
    plotdir = '_plots'
    anim = animation_tools.make_anim(plotdir, fname_pattern='fgout*.png',
                figsize=(13,8))
    fname = 'animation.html'
    fname_mp4 = f'particle_animation_{event}.mp4'
    title = 'Test Animation'
    #animation_tools.make_html(anim, fname, title, fps=5)
    animation_tools.make_mp4(anim, fname_mp4, fps=5)

if 1:
    # better way to make animation:
    fgout = fgout0
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
        t_str = P.timeformat(t)
        title_text = title('Water depth at %s after quake' % t_str)

        points_water, = P.plot_debris(fgout.t, debris_paths, dbnos_water, ax,
                            color='b', marker='s',markersize=0.5)
        points_land, = P.plot_debris(fgout.t, debris_paths, dbnos_land, ax,
                            color='r', marker='s',markersize=0.5)



    def update(fgframe):
        #fgframe = framenos[num]
        fgout = P.read_fgout_frame(fgno, fgframe, plotdata)
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

        t_str = P.timeformat(t)
        title_text.set_text('Water depth at %s after quake' % t_str)
        #return im,points_water,points_land,title_text

    print('Making anim...')
    anim = animation.FuncAnimation(fig, update,
                          frames=fgframes, interval=200, blit=False)

    fname_mp4 = f'landwater_tracers_{event}.mp4'
    fps = 5
    print('Making mp4...')
    animation_tools.make_mp4(anim, fname_mp4, fps)
