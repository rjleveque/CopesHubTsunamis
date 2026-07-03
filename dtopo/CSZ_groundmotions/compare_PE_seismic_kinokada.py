from pylab import *
from clawpack.geoclaw import topotools, dtopotools
from clawpack.clawutil.util import fullpath_import

CHTtools = fullpath_import('/Users/rjl/git/CHTuser/src/CHTuser/CHTtools.py')

coast = load('CSZ_coast.npy')
etopo = topotools.read_netcdf('etopo22_30sec', extent=[-129,-123,40,50],
                              coarsen = 2)
etopo_fcn = etopo.make_function()

datadir1 = 'dtopo30sec/dtopofiles'
label1 = 'seismic, with subevents'
event1 = lambda event: event

if 0:
    datadir2 = 'dtopo30sec_nosubevents_kinokada/dtopofiles'
    label2 = 'Okada, no subevents'
    event2 = lambda event: event

if 1:
    datadir2 = 'okada/dtopofiles_instant'
    label2 = 'Okada, with subevents'
    event2 = lambda event: CHTtools.longname(event[:5]) + '_okada_instant'

event = 'BL10S_instant'

dtopo1 = dtopotools.DTopography(f'{datadir1}/{event1(event)}.dtt3')
dtopo2 = dtopotools.DTopography(f'{datadir2}/{event2(event)}.dtt3')

#assert abs(dtopo1.x - dtopo2.x).max() < 1e-12, 'dtopo x arrays do not agree'
#assert abs(dtopo1.y - dtopo2.y).max() < 1e-12, 'dtopo y arrays do not agree'

etopo_dtopo1 = etopo_fcn(dtopo1.X, dtopo1.Y)
etopo_dtopo2 = etopo_fcn(dtopo2.X, dtopo2.Y)

figure()
j = where(dtopo1.y < 47)[0].max()
plot(dtopo1.x, dtopo1.dZ[-1,j,:], 'b', label='version 1')
plot(dtopo2.x, dtopo2.dZ[-1,j,:], 'r', label='version 2')

fig,axs = subplots(1,2,figsize=(11,8))
fig.suptitle(event, fontsize=20)
dtopo1.plot_dZ_colors(10,axes=axs[0], dZ_interval=1e6, cmax_dZ=10)
axs[0].set_title(label1)
axs[0].plot(coast[:,0], coast[:,1], 'g', linewidth=0.7)
axs[0].set_xlim(-129, -123)
axs[0].set_ylim(40,50)

dtopo2.plot_dZ_colors(10,axes=axs[1], dZ_interval=1e6, cmax_dZ=10)
axs[1].set_title(label2)
axs[1].plot(coast[:,0], coast[:,1], 'g', linewidth=0.7)
axs[1].set_xlim(-129, -123)
axs[1].set_ylim(40,50)

dx1 = dtopo1.x[1] - dtopo1.x[0]
dy1 = dtopo1.y[1] - dtopo1.y[0]
eta1_offshore = where(etopo_dtopo < 0, dtopo1.dZ[-1,:,:], 0.)
PE1 = 0.5 * (9.81 * 1025 * eta1_offshore**2).sum() * dx1*dy1
print(f'dtopo1 has total potential energy {PE1:.4e} Joules')

dx2 = dtopo2.x[1] - dtopo2.x[0]
dy2 = dtopo2.y[1] - dtopo2.y[0]
eta2_offshore = where(etopo_dtopo < 0, dtopo2.dZ[-1,:,:], 0.)
PE2 = 0.5 * (9.81 * 1025 * eta2_offshore**2).sum() * dx2*dy2
print(f'dtopo2 has total potential energy {PE2:.4e} Joules')

# all events:

depths = ['D','M','S']

# buried_locking events:
all_events = [f'BL10{depth}' for depth in depths] \
           + [f'BL13{depth}' for depth in depths] \
           + [f'BL16{depth}' for depth in depths] \

# add random events:
all_events += [e.replace('L','R') for e in all_events]

# add ft events:
#all_events += [e.replace('B','F') for e in all_events]

all_events.sort()
all_events = [e + '_instant' for e in all_events]

PE1_vals = []
PE2_vals = []
for event in all_events:
    dtopo1 = dtopotools.DTopography(f'{datadir1}/{event1(event)}.dtt3')
    dtopo2 = dtopotools.DTopography(f'{datadir2}/{event2(event)}.dtt3')
    eta1_offshore = where(etopo_dtopo < 0, dtopo1.dZ[-1,:,:], 0.)
    eta2_offshore = where(etopo_dtopo < 0, dtopo2.dZ[-1,:,:], 0.)
    PE1 = 0.5 * (9.81 * 1025 * eta1_offshore**2).sum() * dx1*dy1
    PE2 = 0.5 * (9.81 * 1025 * eta2_offshore**2).sum() * dx2*dy2
    PE1_vals.append(PE1)
    PE2_vals.append(PE2)

PE1_vals = array(PE1_vals)/1e3
PE2_vals = array(PE2_vals)/1e3

ind_PE = argsort(PE2_vals)
print('Largest events based on PE2:')
for i in range(len(ind_PE)-1,-1,-1):
    k = ind_PE[i]
    print(f'{all_events[k]}:  {PE1_vals[k]:7.3f}  {PE2_vals[k]:7.3f}')

figure()
PEmax = 800
plot(PE1_vals[:18], PE2_vals[:18], 'bo', markersize=4)
plot(PE1_vals[18:], PE2_vals[18:], 'ro', markersize=4)
plot([0,PEmax], [0, PEmax], 'k-')
axis('scaled')
xlim(0,PEmax)
ylim(0,PEmax)
grid(True)
xlabel(f'Potential Energy of {label1} (kJ)')
ylabel(f'Potential Energy of {label2} (kJ)')
title('Scatter plot of Potential Energy (in kiloJoules)\n' \
    + 'for buried (blue) and frontal thrust (red) events')

for i in range(len(ind_PE)-1, len(ind_PE)-6, -1):
    k = ind_PE[i]
    text(PE1_vals[k], PE2_vals[k]+10, all_events[k][:5],
         ha='center',fontsize=9)
