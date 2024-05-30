"""
Adapted from CMH project to show ship tracks
"""

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend


from pylab import *
import os,sys
from matplotlib import colors
import matplotlib.animation as animation
from clawpack.visclaw import animation_tools, plottools, geoplot
from clawpack.geoclaw import topotools


try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Set CHT enviornment variable to repository top")
    
sys.path.insert(0,os.path.join(CHT,'common_code'))
import fgout_particles as P


if 1:
    from clawpack.geoclaw import fgout_tools
    graphics_dir = os.path.abspath('../../../graphics')
else:
    # local versions for self-contained directory:
    import fgout_tools
    graphics_dir = './'

#topo = topotools.Topography()
#topo.read('../BM3topo.tt3',3)

#graphics_dir = '/Users/rjl/git/CMH_2021/graphics/'

fgno = 3

# List of frames to use for making debris paths and animation:
#fgframes = range(100,110)
#fgframes = range(1,331)
#fgframes = [100,200]
fgframes = range(1,482,2)

fgout_extent = [-123.47,-123.36,48.12,48.15]
plot_extent = fgout_extent

#fgout_extent = [-123.85, -123.81, 46.185, 46.196]
#plot_extent = [-123.84, -123.81, 46.185, 46.196]

#bgimage = [imread('fgout1GE.jpg'), fgout_extent]
#bgimage = [imread(os.path.join(graphics_dir,'fgout1STfootprint.png')),
#           plot_extent]
bgimage = None
name = 'ships'

outdir = '_output'
format = 'binary32'  # format of fgout grid output


color = 'k'
linewidth = 3


# Create ClawPlotData object used for reading in fgout frames:
#plotdata = P.make_plotdata(fgno, outdir, format)

fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format) 

fgout100 = fgout_grid.read_frame(100)
Bfine = fgout100.B

# Deterime time t0 of first fgout frame, to initialize particles
frameno0 = fgframes[0]
fgout0 = fgout_grid.read_frame(frameno0)
t0 = fgout0.t


def timeformat(t):
    """
    Convert t in seconds to string in format hh:mm:ss
    """
    from numpy import mod
    hours = int(t/3600.)
    tmin = mod(t,3600.)
    min = int(tmin/60.)
    sec = int(mod(tmin,60.))
    timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
    return timestr
    
    
# Initialize debris_paths dictionary and set
# initial debris particle locations (x0,y0) at time t0.
# Require a list of dbnos and each array 
#     debris_paths[dbno]
# in the dictionary is a 2d array with a single row [t0, x0, y0] to start.

debris_paths = {}
grounding_depth = {}
drag_factor = {}

dbnos = []
dbnosA = []
dbnosB = []
# set initial velocities to 0 (or may want to interpolate from fgout0 if t0>0?)
u0 = 0.
v0 = 0.

        
grounding_depth_common = 5.
drag_factor_common = None

ship_ends = [[-123.4178, 48.1395, -123.4168, 48.1396],
             [-123.4140, 48.1395, -123.4130, 48.1396]
            ]
for k in range(len(ship_ends)):
    dbno = k+1
    x1,y1,x2,y2 = ship_ends[k]
    db = array([[t0, x1, y1, u0, v0]])
    debris_paths[dbno] = db
    dbnos.append(dbno)
    dbnosA.append(dbno)
    grounding_depth[dbno] = grounding_depth_common
    drag_factor[dbno] = drag_factor_common 
    
    # paired particle:
    dbno = 100+dbno
    db = array([[t0, x2, y2, u0, v0]])
    debris_paths[dbno] = db
    dbnos.append(dbno)
    dbnosB.append(dbno)        
    grounding_depth[dbno] = grounding_depth_common
    drag_factor[dbno] = drag_factor_common
    
nships = len(dbnos)

dbnos_water = []
if 0:
    # still need to fix...
    x1 = -123.849 # 0.01 west of fgout
    x2 = -123.813 # east of lake
    y1 = 46.189
    y2 = 46.196
    dxd = 0.0004
    
    dyd = dxd * cos(47*pi/180)
    xd = arange(x1, x2, dxd)
    yd = arange(y1, y2, dyd)
    mxd = len(xd)
    myd = len(yd)
    print('Will put down %i by %i array of particles' % (mxd,myd))

    for i in range(len(xd)):
        x0 = xd[i]
        for j in range(len(yd)):
            y0 = yd[j]
            dbno = nships + myd*i + j
            db = array([[t0, x0, y0, u0, v0]])
            debris_paths[dbno] = db
            
            if B_fcn_xy(x0,y0) < -1:
                dbnos_water.append(dbno)
                grounding_depth[dbno] = 0.
                drag_factor[dbno] = None
            if 0:  # B_fcn_xy(x0,y0) > -1:
                dbnos_land.append(dbno)
                grounding_depth[dbno] = 2.
                drag_factor[dbno] = 0.1

    # does this work??
    debris_paths_water = P.make_debris_paths(fgout_grid, fgframes, 
                                       debris_paths, dbnos_water,
                                       drag_factor, grounding_depth)
    
print('Created %i initial debris particles' % len(dbnos))


                     
# Compute debris path for each particle by using all the fgout frames
# in the list fgframes (first frame should be frameno0 used to set t0 above):

debris_paths = P.make_debris_paths_pairs(fgout_grid, fgframes,  
                                         debris_paths, dbnosA, dbnosB,
                                         drag_factor=drag_factor, 
                                         grounding_depth=grounding_depth)


def make_dbAB(t, debris_paths, dbnosA, dbnosB):
    xdAB = []
    ydAB = []
    
    for k in range(len(dbnosA)):
        dbnoA = dbnosA[k]
        dbnoB = dbnosB[k]
        dbA = debris_paths[dbnoA]
        dbB = debris_paths[dbnoB]
        try:
            j = where(abs(dbA[:,0]-t) < 1e-6)[0].max()
        except:
            print('Did not find path for dbno=%i at t = %.3f' % (dbno,t))
            j = -1
        if j > -1:
            xA = dbA[j,1]
            yA = dbA[j,2]
            xB = dbB[j,1]
            yB = dbB[j,2]
            xdAB = xdAB + [xA,xB,nan]
            ydAB = ydAB + [yA,yB,nan]
            #import pdb; pdb.set_trace()
    xdAB = array(xdAB)
    ydAB = array(ydAB)
    return xdAB,ydAB



# First initialize plot with data from initial frame,
# do this in a way that returns an object for each plot attribute that
# will need to be changed in subsequent frames.
# In tis case, the color image of the water depth, plots of particles, and 
# title (which includes the time) will change.
# The background image, colorbar, etc. do not change.

fgout = fgout0
ylat = fgout.Y.mean()  # for aspect ratio of plots

fig,ax = subplots(figsize=(12,7))
if bgimage:
    ax.imshow(bgimage[0],extent=bgimage[1])
    
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])


s_units = 'knots'
if s_units == 'knots':
    s = fgout.s * 1.9438  # convert m/s to knots
    bounds_speed = np.array([1e-3,2,4,6,8,10])  # knots
    #bounds_speed = np.array([1e-3,1,2,3,6,9])  # knots
elif s_units == 'm/s':
    s = fgout.s
    bounds_speed = np.array([1e-3,0.05,0.1,0.2,0.3,0.4])  # m/s
elif s_units == 'vorticity':
    dx = (fgout.x[1] - fgout.x[0]) * 111e3 * cos(ylat*pi/180)
    dy = (fgout.y[1] - fgout.y[0]) * 111e3
    #print('dx = %.2fm, dy = %.2fm' % (dx,dy))
    dudy = zeros(fgout.u.shape)
    dvdx = zeros(fgout.u.shape)
    dudy[:,:-1] = diff(fgout.u, axis=1) / dy
    dvdx[:-1,:] = diff(fgout.v, axis=0) / dx
    s = dvdx - dudy
    omega_max = abs(s).max()
    print('+++ omega_max = %.1f' % omega_max)
    bounds_speed = np.array([-3,-2,-1,0,1,2])
    
#a = 0.4 # on top of image
a = 0.5
cmap_speed = mpl.colors.ListedColormap([
                [.3,.3,1,a],[0,0,1,a], [1,.8,.8,a],[1,.6,.6,a],
                [1,0,0,a]])

# Set color for value exceeding top of range to purple:
cmap_speed.set_over(color=[1,0,1,a])

if bgimage:
    # Set color to transparent where s < 1e-3:
    cmap_speed.set_under(color=[1,1,1,0])
else:
    # Set color to white where s < 1e-3:
    cmap_speed.set_under(color=[1,1,1,a])

norm_speed = colors.BoundaryNorm(bounds_speed, cmap_speed.N)

#contourf(fgout.X, fgout.Y, fgout.B, levels=[0.04, 0.2], colors=[[1,1,.7]])
#Bdry = where(fgout.h < 0.001, fgout.B, nan)
#pcolormesh(fgout.X, fgout.Y, Bdry, cmap=geoplot.land1_colormap)
im = imshow(flipud(s.T), extent=fgout_extent,
            cmap=cmap_speed, norm=norm_speed)
cb = colorbar(im, extend='max', shrink=0.7)
cb.set_label(s_units)
#contour(fgout.X, fgout.Y, Bfine, [0], colors='g', linewidths=1)
#contour(fgout.X, fgout.Y, Bfine, [0.02,0.04,0.08], colors='g', linewidths=0.5)


ax.set_aspect(1./cos(ylat*pi/180.))
ticklabel_format(useOffset=False)
xticks(rotation=20)
ax.set_xlim(plot_extent[:2])
ax.set_ylim(plot_extent[2:])

t = fgout.t
t_str = timeformat(t)
title_text = title('Water speed at %s after quake' % t_str)
#title_text = title('Vorticity at %s after quake' % t_str)


xdAB,ydAB = make_dbAB(t, debris_paths, dbnosA, dbnosB)
pairs, = ax.plot(xdAB, ydAB, color=color, linestyle='-', linewidth=linewidth)
print('+++ pairs = ',pairs)

# The function update below should have arguments num (for the frame number)
# plus things listed here in fargs.

fargs = (im,pairs,title_text)

# fargs should be initialized above and are the plot Artist objects 
# whose data change from one frame to the next.

if 0:
    # Plot full debris paths on the previous plot:
    #for k in range(len(dbnos)):
    for k in [6,7]:
        debris_path = debris_paths[dbnos[k]]
        ax.plot(debris_path[:,1], debris_path[:,2], color='k', linewidth=0.5)

def update(num, im, pairs, title_text):
    
    fgframe = fgframes[num]
    # note: uses fgframes to specify fgout frames to use
    
    # Read fgout data for this frame:
    #fgout = P.read_fgout_frame(fgno, fgframe, plotdata)
    fgout = fgout_grid.read_frame(fgframe)
    
    # Reset the plot objects that need to change from previous frame:

    # title:
    t = fgout.t        
    t_str = timeformat(t)
    title_text = title('Water speed at %s after quake' % t_str)
    #title_text = title('Vorticity at %s after quake' % t_str)

    
    # color image:
    if s_units == 'knots':
        s = fgout.s * 1.9438  # convert m/s to knots
    elif s_units == 'm/s':
        s = fgout.s
    elif s_units == 'vorticity':
        dx = (fgout.x[1] - fgout.x[0]) * 111e3 * cos(ylat*pi/180)
        dy = (fgout.y[1] - fgout.y[0]) * 111e3
        #print('dx = %.2fm, dy = %.2fm' % (dx,dy))
        dudy = zeros(fgout.u.shape)
        dvdx = zeros(fgout.u.shape)
        dudy[:,:-1] = diff(fgout.u, axis=1) / dy
        dvdx[:-1,:] = diff(fgout.v, axis=0) / dx
        s = dvdx - dudy
        omega_max = abs(s).max()
        print('+++ omega_max = %.1f' % omega_max)
    
    im.set_data(flipud(s.T))

    # particle locations:
    
    xdAB,ydAB = make_dbAB(t, debris_paths, dbnosA, dbnosB)
    pairs.set_data(xdAB, ydAB)        
        
    # must now return all the objects listed in fargs:
    return im,pairs,title_text

print('Making anim...')
anim = animation.FuncAnimation(fig, update,
                               frames=len(fgframes), 
                               fargs=fargs,
                               interval=200, blit=True)

fname_mp4 = '%s_pairs.mp4' % name
fps = 5
print('Making mp4...')
animation_tools.make_mp4(anim, fname_mp4, fps)
