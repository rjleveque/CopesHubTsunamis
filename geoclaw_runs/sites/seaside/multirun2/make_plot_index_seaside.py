"""
Code to make plot indices for  multiple GeoClaw runs,
with a different dtopo file for each.
"""

import os,sys

location = 'Seaside'

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.abspath('..'))
#import plot_gauges, process_fgmax


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

print('geoclaw_plots = ',geoclaw_plots)

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

events.sort()

plotdirs = ['%s/_plots_%s' % (geoclaw_plots, event) for event in events]


def make_html_index(plotdir,event):
    html_fname = os.path.join(plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write('\n<ul>\n<li><a href="fgmax">fgmax plots</a>\n')
        f.write('<li><a href="gauges">gauge plots</a>\n</ul>\n')
    print('Created ',html_fname)

top_index_fname = os.path.join(geoclaw_plots,'index.html')
with open(top_index_fname, 'w') as top_index:
    
    top_index.write('<html>\n<h1>Plots for %s</h1>\n<ul>\n' % location)

    for k in range(len(plotdirs)):
        plotdir = plotdirs[k]
        event = events[k]
        run_name = '%s_%s' % (location,event)

        make_html_index(plotdir,event)
        
        # relative paths:
        plotdir = '_plots_%s' % event
        gauges_plotdir = '_plots_%s/gauges' % event
        fgmax_plotdir = '_plots_%s/fgmax' % event

        top_index.write('<h2>%s</h2>\n ' % event)

        # images to appear on index page:

        if 1:
            # fgmax plots
            top_index.write('<img src=%s/%s_h_onshore.png width=28%%>\n' \
                    % (fgmax_plotdir,run_name))
            top_index.write('<img src=%s/%s_transects.png width=38%%>\n' \
                    % (fgmax_plotdir,run_name))
            top_index.write('<img src=%s/%s_speed.png width=28%%>\n<p>\n' \
                    % (fgmax_plotdir,run_name))                
    
            top_index.write('<img src=%s/%s_Gauge01033.png width=30%%>\n' \
                    % (gauges_plotdir,run_name))
            top_index.write('<img src=%s/%s_Gauge01045.png width=30%%>\n' \
                    % (gauges_plotdir,run_name))
            top_index.write('<img src=%s/%s_Gauge01047.png width=30%%>\n<p>\n' \
                    % (gauges_plotdir,run_name))
                
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
            top_index.write('<li> <a href="%s">gauges</a> On 12th Avenue:\n' \
                    % gauges_plotdir)
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_Gauge01033.png">1033 (offshore)</a>\n' \
                    % (gauges_plotdir,run_name))    
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_Gauge01013.png">1013 (in Necamicum River)</a>\n' \
                    % (gauges_plotdir,run_name)) 
            top_index.write('&nbsp;&nbsp; <a href="%s/%s_Gauge01017.png">1017 (in Neawanna River)</a>\n' \
                    % (gauges_plotdir,run_name))                

        top_index.write('</ul>\n')
        
print('Created', top_index_fname)
