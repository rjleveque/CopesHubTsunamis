"""
Code to make plot indices for  multiple GeoClaw runs,
with a different dtopo file for each.
"""

import os,sys

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.abspath('..'))
from make_plot_index import make_index  # eventually move to common_code?


location = 'Newport'
show_gaugenos = [1007, 1022, 1041]

# set geoclaw_plots to point to directory containing plots for each event:

# location for big files:
this_dir = os.getcwd()
# Randy's laptop:
scratch_dir = this_dir.replace('git/CopesHubTsunamis/geoclaw_runs', \
                               'scratch/CHT_runs')
# for hyak:
scratch_dir = scratch_dir.replace('/mmfs1/home', '/gscratch/tsunami')

# where to find output for all the runs:
# (in a subdirectory of runs_dir named geoclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)
#runs_dir = os.path.abspath('.')

# top level of plots directory:
geoclaw_plots = os.path.join(runs_dir, 'geoclaw_plots')


geoclaw_plots = './geoclaw_plots'  # for local version
print('geoclaw_plots = ',geoclaw_plots)


# Select set of events to show plots for:

if 0:

    #CoPes Hub ground motions:

    all_models = \
        ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
         'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']

    if 1:
        models = all_models
        #models = all_models[:3]
        events = ['%s-deep' % model for model in models] \
               + ['%s-middle' % model for model in models] \
               + ['%s-shallow' % model for model in models]

    instant = True
    if instant:
        events = [e+'_instant' for e in events]
        rupture_type = 'Instant'

    if 0:
        #events = ['buried-random-str10-middle','buried-random-str10-shallow']
        events = ['buried-locking-str10-deep']
        rupture_type = 'Kinematic'

    events.sort()


# for Tshirts:

sizes = ['SM','M','L','XL','XXL']

all_events = []
for size in sizes:
    for M in [1,2,3]:
        all_events.append('CSZ_%s%s_noext' % (size,M))

events = []
for e in all_events:
    if e[4] in ['M','L']: events.append(e)

#events = ['CSZ_M1_noext', 'CSZ_L1_noext']



plotdirs = ['%s/_plots_%s' % (geoclaw_plots, event) for event in events]

cmd = 'cp ../%sGauges.* %s' % (location,geoclaw_plots)
print(cmd)
os.system(cmd)


def make_html_index(plotdir,event,gauges_plotdir,fgmax_plotdir):
    """
    Make index for a single event.
    This function is called in loop below.
    """
    html_fname = os.path.join(plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write('\n<ul>\n<li><a href="fgmax">fgmax plots</a>\n')
        f.write('<li><a href="gauges">gauge plots</a>\n</ul>\n')
    print('Created ',html_fname)


# top_index show select plots from each event with links to plot index for each
top_index_fname = os.path.join(geoclaw_plots,'index.html')

with open(top_index_fname, 'w') as top_index:

    #top_index.write('<html>\n<h1>Plots for %s</h1>\n<ul>\n' % location)
    top_index.write('<html>\n<h1>Inundation plots for %s</h1>\n' % location)
    top_index.write(
"""Computed from T-shirt size events
using the GeoClaw tsunami model.\n
""")
    top_index.write(
"""<h2> Gauge locations:</h2>
Download <a href="%sGauges.kml">%sGauges.kml</a> and open in
Google Earth.
<p>
<img src=%sGauges.jpg width=50%%>
<p>
&nbsp;
<p>
""" % (location,location,location))

    for k in range(len(plotdirs)):
        plotdir = plotdirs[k]
        event = events[k]
        run_name = '%s_%s' % (location,event)

        gauges_plotdir = gauges'
        fgmax_plotdir = 'fgmax/fgmax1' # Yaquina Bay

        make_html_index(plotdir,event,gauges_plotdir,fgmax_plotdir)

        # relative paths:
        plotdir = '_plots_%s' % event

        gauges_plotdir = '%s/%s' % (plotdir, gauges_plotdir)
        fgmax_plotdir = '%s/%s' % (plotdir, fgmax_plotdir)

        top_index.write('<h2>%s</h2>\n ' % event)

        # images to appear on index page:

        # fgmax plots
        top_index.write('<img src=%s/%s_h_onshore.png width=28%%>\n' \
                % (fgmax_plotdir,run_name))
        top_index.write('<img src=%s/%s_transects.png width=38%%>\n' \
                % (fgmax_plotdir,run_name))
        top_index.write('<img src=%s/%s_speed.png width=28%%>\n<p>\n' \
                % (fgmax_plotdir,run_name))

        # select gauges:
        for gaugeno in show_gaugenos:
            top_index.write('<img src=%s/%s_Gauge%s.png width=30%%>\n' \
                    % (gauges_plotdir,run_name,str(gaugeno).zfill(5)))


        # Other links to appear in index:

        top_index.write('<ul>\n<li> <a href="%s/index.html">all plots</a>\n' \
                % plotdir)
        top_index.write('<li>  animation: <a href="%s/%s_animation.mp4">mp4</a>\n' \
                % (plotdir,run_name))

        if 1:

            top_index.write('<li>  fgmax: <a href="%s/%s_h_onshore.png">max h onshore</a>\n' \
                    % (fgmax_plotdir,run_name))
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_speed.png">max speed</a>\n' \
                    % (fgmax_plotdir,run_name))
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_transects.png">transects</a>\n' \
                    % (fgmax_plotdir,run_name))
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_fgmax.kmz">kmz file</a>\n' \
                    % (fgmax_plotdir,run_name))

            top_index.write('<li> <a href="%s">gauges</a>Gauges:\n' \
                    % gauges_plotdir)


        top_index.write('</ul>\n')

print('Created', top_index_fname)
