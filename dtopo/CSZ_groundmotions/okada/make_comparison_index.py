"""
Code to make plot indices for  multiple GeoClaw runs,
with a different dtopo file for each.
"""

import os,sys

#location = 'offshore'

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

#sys.path.insert(0, os.path.abspath('..'))
#import plot_gauges, process_fgmax



all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']

if 1:
    models = all_models
    #models = all_models[:3]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

if 0:
    #events = ['buried-random-str10-middle','buried-random-str10-shallow']
    events = ['buried-locking-str10-deep']


plotdir = 'plots_compare_seismic_okada_instant'


top_index_fname = os.path.join(plotdir,'index.html')
with open(top_index_fname, 'w') as top_index:

    top_index.write('<html>\n<h1>Comparisons of Seismic and Okada (static - final vertical displacements)</h1>\n<ul>\n')

    for event in events:
        cplot = 'compare_okada_seismic_instant_%s.png' % event
        tplot = 'transects_okada_seismic_instant_%s.png' % event

        top_index.write('<b>%s</b><p>\n' % event)
        top_index.write('<img src=%s width=600>\n' % cplot)
        top_index.write('<img src=%s width=500>\n<p>\n' % tplot)


print('Created', top_index_fname)
