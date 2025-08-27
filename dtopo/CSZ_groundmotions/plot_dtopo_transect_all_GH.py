
"""
Plot transects of dtopo deformations at final time
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools,dtopotools

# transect latitude:
#y0 = 46.
y0_list = [46.95]

all_models = []

if 0:
    all_models = all_models + \
        ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
         'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']
    name_png = 'transects_GH_buried.png'

if 1:
    all_models = all_models + \
        ['ft-locking-mur13', 'ft-locking-skl16', 'ft-locking-str10',
         'ft-random-mur13',  'ft-random-skl16',  'ft-random-str10']
    name_png = 'transects_GH_ft.png'

if len(all_models) == 12:
    # including both buried and ft:
    name_png = 'transects_GH_all.png'

models = all_models
#models = all_models[:3]
events = ['%s-deep' % model for model in models] \
       + ['%s-middle' % model for model in models] \
       + ['%s-shallow' % model for model in models]

events.sort()

for y0 in y0_list:

    figure(205,figsize=(10,7))
    clf()
    #linecolors = ['r','b','g']
    #labels =

    for k,event in enumerate(events):
        print('Rupture name: ',event)
        rupture = event

        # Read dtopo file used for GeoClaw:

        fname_dtopo = 'dtopofiles_FrontalThrust/' + event + '.dtt3'
        print('Reading ',fname_dtopo)
        dtopo = dtopotools.DTopography(fname_dtopo, 3)

        j = where(dtopo.y<y0)[0].max()
        #plot(dtopo.x,dtopo.dZ[-1,j,:],linecolors[k],label=labels[k])
        plot(dtopo.x,dtopo.dZ[-1,j,:],label=event)
        title('Grays Harbor\nTransect of final vertical displacement at y = %.2f' \
                % y0,fontsize=15)
    grid(True)
    #xlim(-127,-122)
    xlim(-124.2, -123.7)
    #xlim(-124.2, -123.3)
    ylim(-4,2)

    if 1:
        events = []
        depths = ['L1','L3','M1']
        rupture = 'CSZ_'
        for depth in depths:
            events.append('%s%s_noext' % (rupture,depth))

        linecolors = ['k','m','c']
        labels = depths

        for k,event in enumerate(events):
            print('Rupture name: ',event)

            fname_dtopo = '../CSZ_Tshirts/' + event + '.tt3'
            print('Reading ',fname_dtopo)
            dtopo = dtopotools.DTopography(fname_dtopo, 3)

            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolors[k],
                 linestyle='--', label=labels[k])

    if 0:
        events = []
        labels = ['Region2','Region3']
        for label in labels:
            events.append('ASCE_SIFT_dtopo_%s' % label)

        #linecolors = ['m','c','g']
        linecolors = ['k','m','g']

        for k,event in enumerate(events):
            print('Rupture name: ',event)

            fname_dtopo = '/Users/rjl/B/dtopo/ASCE_SIFT/' + event + '.tt3'
            print('Reading ',fname_dtopo)
            dtopo = dtopotools.DTopography(fname_dtopo, 3)

            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolors[k],
                 linestyle='--', label=labels[k])

    if 0:
        events = []
        labels = ['S-A-Whole','B-Whole','f1p4as1']
        for label in labels[:2]:
            events.append('%s-result' % label)
        events.append('f1p4as1-interpolated-new')

        #linecolors = ['m','c','g']
        linecolors = ['k','m','c']

        for k,event in enumerate(events):
            print('Rupture name: ',event)

            fname_dtopo = '/Users/rjl/B/dtopo/CSZ/Wang/dtopofiles/' + event + '.tt3'
            print('Reading ',fname_dtopo)
            dtopo = dtopotools.DTopography(fname_dtopo, 3)

            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolors[k],
                 linestyle='--', label=labels[k])

    #legend(loc='upper right',fontsize=12,framealpha=1)

    if 1:
        fname = 'dtopo_transect_GH_ft.png'
        savefig(fname)
        print('Created ',fname)
