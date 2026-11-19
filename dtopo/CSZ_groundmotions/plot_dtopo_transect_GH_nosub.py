
"""
Plot transects of dtopo deformations at final time

Comparing original seismic with new nosub_kinokada version
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools,dtopotools

# transect latitude:
y0 = 46.95

plotdir = 'plots_GH_SeismicVsNosub_46-95'


try:
    import CHTtools
    all_events = CHTtools.all_events()  # returns the list of 36 events
except:
    # if CHTtools isn't found, construct the list of 36 events here:

    depths = ['D','M','S']

    # buried_locking events:
    all_events = [f'BL10{depth}' for depth in depths] \
               + [f'BL13{depth}' for depth in depths] \
               + [f'BL16{depth}' for depth in depths] \

    all_events += [e.replace('L','R') for e in all_events]  # add random events
    all_events += [e.replace('B','F') for e in all_events]  # add ft events

    all_events.sort()

events = all_events


add_Tevents = True

if add_Tevents:
    Tevents = ['M1','L3']

    Tdtopos = {}
    for k,Tevent in enumerate(Tevents):
        print('Rupture name: ',Tevent)

        fname_dtopo = f'../CSZ_Tshirts/CSZ_{Tevent}_noext.tt3'
        print('Reading ',fname_dtopo)
        Tdtopos[Tevent] = dtopotools.DTopography(fname_dtopo, 3)

    Tj = where(Tdtopos[Tevent].y<y0)[0].max()

for event in events:

    figure(205,figsize=(12,7))
    clf()

    fname_dtopo = 'dtopo30sec_seismic/dtopofiles/' + event + '.dtt3'
    print('Reading ',fname_dtopo)
    dtopo = dtopotools.DTopography(fname_dtopo, 3)

    j = where(dtopo.y<y0)[0].max()
    plot(dtopo.x,dtopo.dZ[-1,j,:],'r',linewidth=2, label='Original seismic')

    fname_dtopo = 'dtopo30sec_nosubevents_kinokada/dtopofiles_nc/' + event + '.nc'
    print('Reading ',fname_dtopo)
    dtopo = dtopotools.DTopography(fname_dtopo, 4)

    j = where(dtopo.y<y0)[0].max()
    plot(dtopo.x,dtopo.dZ[-1,j,:],'b',linewidth=2, label='No Subevents, KinOkada')


    title(f'Event {event}\n' \
          + f'Transect of final vertical displacement at latitude {y0:.2f}', 
            fontsize=15)
    
    grid(True)
    #xlim(-127,-122)
    #xlim(-124.2, -123.7)
    #xlim(-124.6, -123.7)  # out to 100m isobath
    xlim(-127, -123) 
    dz_min = -4
    dz_max = 8
    ylim(dz_min, dz_max)
    xlabel('Longitude', fontsize=12)
    ylabel('Final vertical deformation (m)', fontsize=12)


    if add_Tevents:
        linecolors = ['k','m','c']
        for k,Tevent in enumerate(Tevents):
            plot(Tdtopos[Tevent].x, Tdtopos[Tevent].dZ[-1,Tj,:],
                 color=linecolors[k],
                 linestyle='--', label=f'T-shirt {Tevent}')


    legend(loc='upper right',fontsize=12,framealpha=1)

    fill_between([-124.18,-124.13], [dz_min,dz_min], [dz_max,dz_max],
                 color=[.9,1,.9])
    text(-124.155, dz_min+0.1, 'Ocean\nShores', color='g', ha='center',
        fontsize=10)

    fill_between([-123.88,-123.80], [dz_min,dz_min], [dz_max,dz_max],
                 color=[.9,1,.9])
    text(-123.84, dz_min+0.1, 'Aberdeen', color='g', ha='center',
        fontsize=10)

    if 1:
        fname = f'{plotdir}/dtopo_transect_GH_SeismicVsNosub_{event}.png'
        savefig(fname)
        print('Created ',fname)
