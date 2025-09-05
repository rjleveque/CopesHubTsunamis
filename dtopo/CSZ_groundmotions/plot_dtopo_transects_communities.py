
"""
Plot transects of dtopo deformations at final time

Refactored to loop over all rupture models, and T-shirt reference events,
reading in each dtopo file only once.
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools,dtopotools

# transect latitude:
#y0_list = [46.95]

transects = [
    (47.91, [(-124.64, 'La Push')]),
    (46.98, [(-124.16,'Ocean Shores'),(-123.86, 'Aberdeen')]),
    (46.0, [(-123.93, 'Seaside')]),
    (44.62, [(-124.07, 'Newport')]),
    (40.80, [(-124.24, 'Eureka')]),
]

plotdir = 'seismic_transect_plots'
os.system('mkdir -p %s' % plotdir)


# reference events to appear on all plots:

events = ['CSZ_M1', 'CSZ_M2']
ref_dtopos = []
linecolors = ['k','m','c']
for ke,event in enumerate(events):
    fname_dtopo = '../CSZ_Tshirts/' + event + '_noext.tt3'
    print('Reading ',fname_dtopo)
    dtopo = dtopotools.DTopography(fname_dtopo, 3)
    ref_dtopos.append((dtopo,event,linecolors[ke],'--'))

if 1:

    all_models = \
        ['ft-locking-mur13', 'ft-locking-skl16', 'ft-locking-str10',
         'ft-random-mur13',  'ft-random-skl16',  'ft-random-str10']

    FrontalThrust = True

    if not FrontalThrust:
        all_models = [s.replace('ft','buried') for s in all_models]
        dtopo_dir = 'dtopofiles'
    else:
        dtopo_dir = 'dtopofiles_FrontalThrust'

    ruptures = all_models
    ruptures = all_models[:1]

    ruptures.sort()

for kr,rupture in enumerate(ruptures):
    print('Rupture name: ',rupture)

    events = ['%s-deep' % rupture,
              '%s-middle' % rupture,
              '%s-shallow' % rupture]

    linecolors = ['k','m','c']  # one for each depth

    # loop over transects to set up plots:
    for kt,transect in enumerate(transects):
        y0 = transect[0]
        figure(205+kt,figsize=(13,6))
        clf()
        grid(True)
        dzmin = -6
        dzmax = 16
        xlim(-127,-122)
        ylim(dzmin, dzmax)
        xlabel('Longitude')
        ylabel('Vertical deformation dz (m)')

        title('%s\nTransect of final vertical displacement at y = %.2f' \
                % (rupture,y0),fontsize=15)

        # add dashed lines and labels for coastal location(s):
        for k,clon in enumerate(transect[1]):
            plot([clon[0],clon[0]], [dzmin,dzmax], 'g--',linewidth=1.5)
            text(clon[0], 10.5-k*5, 'Coast at \n%s' % clon[1], ha='center',
                 color='g',fontsize=12,backgroundcolor='w')
        # plot any reference dtopo's:
        for dtopo,label,linecolor,linestyle in ref_dtopos:
            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolor,
                 linestyle=linestyle,label=label)


    # loop over 3 events for this rupture:
    linecolors = ['r','b','g']
    for ke,event in enumerate(events):
        # Read dtopo file used for GeoClaw
        # (use instant version since we only plot final dz, loads faster)
        fname_dtopo = '%s/%s_instant.dtt3' % (dtopo_dir, event)
        print('Reading ',fname_dtopo)
        dtopo = dtopotools.DTopography(fname_dtopo, 3)

        # loop over trasect plots and add curve for this event:
        for kt,transect in enumerate(transects):
            y0 = transect[0]
            figure(205+kt)
            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],linecolors[ke],label=event)

    # loop over transects to add legends, finalize, and save:
    for kt,transect in enumerate(transects):
        y0 = transect[0]
        figure(205+kt)
        legend(loc='upper right',fontsize=10,framealpha=1)

        if 1:
            fname = '%s/dtopo_transects_%s_y%.0f.png' % (plotdir,rupture,100*y0)
            savefig(fname)
            print('Created ',fname)
