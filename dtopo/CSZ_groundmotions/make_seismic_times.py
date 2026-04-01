from pylab import *
import os,sys
from clawpack.geoclaw import dtopotools,util
from clawpack.visclaw import animation_tools, colormaps
from importlib import reload
import OkadaFromSubsampledSlipNosub
reload(OkadaFromSubsampledSlipNosub)


from clawpack.clawutil.util import fullpath_import
CHTuser = os.environ['CHTuser']
CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')

#event = 'BR13D'
event = 'BL10D'


print(f'Loading {event}...')
# just use megathrust even for ft events
fault = OkadaFromSubsampledSlipNosub.load_fault('M',event)

subfaults = [s for s in fault.subfaults if abs(s.slip) > 1e-3]
#rupt_times = array([s.rupture_time for s in fault.subfaults])
latitudes = array([s.corners[0][1] for s in subfaults])
ind = argsort(latitudes)
subfaults = array(subfaults)[ind]
latitudes = latitudes[ind]
slips = array([s.slip for s in subfaults])
rupt_times = array([s.rupture_time for s in subfaults])

fault.rupture_type = 'kinematic'

# Grid for dtopo vertical deformations:
#dx = dy = 30/3600.  # spatial resolution for dtopo file
dx = dy = 2./60.  # for testing
# use same x,y as in dtopo files made from ground motions:
x = arange(-128.5, -122.0001, dx)
y = arange(40,50.0001,dy)
X,Y = meshgrid(x,y,indexing='xy') # indexing used in dtopo

#lat_breaks = array([40, 42, 44, 46, 47.2, 48.3, 50])
lat_breaks = arange(latitudes.min()+0.2, latitudes.max(), 0.2)
subfaults_group = {}
ind_group = {}
ta_first_group = {}
ta_last_group = {}
kgroups = list(range(len(lat_breaks)-1))

# estimates of min, max of seismic wave speeds:
Vmin = 5000
Vmax = 6000

for k in kgroups:
    ind_group[k] = where(logical_and(lat_breaks[k] <= latitudes,
                                     latitudes <= lat_breaks[k+1]))[0]
    #subfaults_group[k] = [s for s in subfaults if \
    #                      (lat_breaks[k] <= s.corners[0][1] < lat_breaks[k+1])]
    subfaults_group[k] = subfaults[ind_group[k]]
    print(f'len(subfaults_group[{k}]) = {len(subfaults_group[k])}')

    ta_first = inf * ones(X.shape)  # accumulate time of first arrival
    ta_last = zeros(X.shape)        # accumulate time of last arrival

    for s in subfaults_group[k][::10]:
        xs, ys, depth = s.corners[0]  # lon,lat,depth of first corner
        xydist = util.haversine(xs, ys, X, Y)  # dist on sphere from fault
        dist = sqrt(xydist**2 + depth**2)  # 3D distance from fault to X,Y point
        ta_first = minimum(ta_first, s.rupture_time + dist/Vmax)
        ta_last = maximum(ta_last, s.rupture_time + s.rise_time + dist/Vmin)

    ta_first_group[k] = ta_first
    ta_last_group[k] = ta_last


coast = load('CSZ_coast.npy')
x_coast = coast[:,0]
y_coast = coast[:,1]


if 0:
    k = 3
    cs = contour(X, Y, ta_first_group[k], colors='k')
    clabel(cs)
    cs = contour(X, Y, ta_last_group[k], colors='b')
    clabel(cs)


figs = []
for t in range(10,700,5):
    fig,ax = subplots(figsize=(6,8))

    plot(x_coast,y_coast,'g',linewidth=0.7)
    axis('scaled')
    axis([-128.5,-122,40,50])

    fault.plot_subfaults(axes=ax,slip_color=True,cmap_slip=colormaps.white_red,
                         slip_time=t,plot_box=False,cmin_slip=0.01,
                         colorbar_shrink=0.7)

    wfrac_groups = {}
    for k in kgroups:

        wfrac = maximum(0, t - ta_first_group[k]) \
                 / (ta_last_group[k] - ta_first_group[k])
        wfrac = minimum(wfrac, 1.)
        wfrac_groups[k] = wfrac
        #contourf(X, Y, wfrac, [.1,.5,.9], colors=[[1,.8,.8],[.8,1,.8]])
        #decay_color = 3*[max(1, t - ta_last_group[k].max())]
        #contour(X, Y, wfrac, [0.33,0.67], colors=[.5,.5,.5], linewidths=0.5)
        contourf(X, Y, wfrac, [0.4,0.5], colors=[[.5,.5,1,.2]])

    if 0:
        for k in kgroups:
            contourf(X, Y, wfrac_groups[k], [0.4,0.6], colors=[[1,.8,.8]])



    figs.append(fig)
    close(fig)

anim = animation_tools.animate_figs(figs, figsize=(6,8))

# make mp4 file:
fname = f'test_seismic_{event}.mp4'
animation_tools.make_mp4(anim, file_name=fname)
print('Created ',fname)



if 0:
    for k in kgroups:
        print(f'k = {k}, ta_first min = {ta_first_group[k].min()}')
        lats = array([s.corners[0][1] for s in subfaults_group[k]])
        print(f'        lat min = {lats.min()}, lat max = {lats.max()}')
        rupt_times = array([s.rupture_time for s in subfaults_group[k]])
        print(f'        rupt_time min = {rupt_times.min()}, rupt_time max = {rupt_times.max()}')
