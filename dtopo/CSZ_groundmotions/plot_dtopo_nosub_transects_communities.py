
"""
Plot transects of dtopo deformations at final time
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

for k,transect in enumerate(transects):
    y0 = transect[0]


    #fname = rupture + '_ASCESIFT_y%s.png' % str(y0).replace('.','-')
    #fname = rupture + '_Wang_y%s.png' % str(y0).replace('.','-')

    figure(305+k,figsize=(13,6))
    clf()

    if 1:
        events = []
        #depths = ['L1','L2','L3']
        depths = ['L1']
        rupture = 'CSZ_'
        for depth in depths:
            events.append('%s%s_noext' % (rupture,depth))

        linecolors = ['k','m','c']
        labels = ['%s%s' % (rupture,dd) for dd in depths]

        for k,event in enumerate(events):
            print('Rupture name: ',event)

            fname_dtopo = '../CSZ_Tshirts/' + event + '.tt3'
            print('Reading ',fname_dtopo)
            dtopo = dtopotools.DTopography(fname_dtopo, 3)

            j = where(dtopo.y<y0)[0].max()
            plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolors[k],
                 linestyle='-.', label=labels[k])

    events = []
    depths = ['deep','middle','shallow']
    rupture = 'buried-locking-str10'
    #rupture = 'buried-locking-mur13'


    rupture = rupture.replace('buried-','buried_')
    events = []
    for depth in depths:
        events.append('%s-%s_NOSUB_SCALED_okada_instant' % (rupture,depth))

    rupture_name = '%s_NOSUB_SCALED_okada' % rupture

    #events = events[:3]
    linecolors = ['r','b','g']
    #labels = depths
    labels = ['%s-%s' % (rupture,depth) for depth in depths]

    for k,event in enumerate(events):
        print('Rupture name: ',event)

        # Read dtopo file used for GeoClaw:

        fname_dtopo = 'audrey_250813_nosub/dtopofiles/' + event + '.dtt3'
        print('Reading ',fname_dtopo)
        dtopo = dtopotools.DTopography(fname_dtopo, 3)

        j = where(dtopo.y<y0)[0].max()
        plot(dtopo.x,dtopo.dZ[-1,j,:],color=linecolors[k],linestyle='-',
             label=labels[k])
        title('%s\nTransect of final vertical displacement at y = %.2f' \
                % (rupture_name,y0),fontsize=15)

    grid(True)
    dzmin = -6
    dzmax = 16
    xlim(-127,-122)
    ylim(dzmin, dzmax)
    xlabel('Longitude')
    ylabel('Vertical deformation dz (m)')

    for kt,clon in enumerate(transect[1]):
        plot([clon[0],clon[0]], [dzmin,dzmax], 'g--',linewidth=1.5)
        text(clon[0], 10.5-kt*5, 'Coast at \n%s' % clon[1], ha='center',
             color='g',fontsize=12,backgroundcolor='w')



    legend(loc='upper right',fontsize=10,framealpha=1)

    if 1:
        plotdir = 'audrey_250813_nosub/transect_plots_nosub'
        os.system('mkdir -p %s' % plotdir)
        fname = '%s/dtopo_transects_%s_y%.0f.png' % (plotdir,rupture_name,100*y0)
        savefig(fname)
        print('Created ',fname)
