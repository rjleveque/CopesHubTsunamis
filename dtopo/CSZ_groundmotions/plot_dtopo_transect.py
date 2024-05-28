
"""
Plot transects of dtopo deformations at final time
"""

from pylab import *
import os,sys,glob
from clawpack.geoclaw import topotools,kmltools,dtopotools

# transect latitude:
#y0 = 46.
y0_list = [45,46,47,48]

for y0 in y0_list:
    events = []
    depths = ['deep','middle','shallow']
    rupture = 'buried-locking-str10'
    for depth in depths:
        events.append('%s-%s_instant' % (rupture,depth))

    fname = rupture + '_ASCESIFT_y%s.png' % str(y0).replace('.','-')

    figure(205,figsize=(10,7))
    clf()
    linecolors = ['r','b','g']
    labels = depths

    for k,event in enumerate(events):
        print('Rupture name: ',event)

        # Save dtopo file for GeoClaw:

        fname_dtopo = 'dtopofiles/' + event + '.dtt3'
        print('Reading ',fname_dtopo)
        dtopo = dtopotools.DTopography(fname_dtopo, 3)

        j = where(dtopo.y<y0)[0].max()
        plot(dtopo.x,dtopo.dZ[-1,j,:],linecolors[k],label=labels[k])
        title('%s\nTransect of final vertical displacement at y = %.1f' \
                % (rupture,y0),fontsize=15)
    grid(True)
    xlim(-127,-122)

    if 0:
        events = []
        depths = ['L1','L2','L3']
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

    if 1:
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

    legend(loc='upper right',fontsize=12,framealpha=1)

    if 1:
        savefig(fname)
        print('Created ',fname)
