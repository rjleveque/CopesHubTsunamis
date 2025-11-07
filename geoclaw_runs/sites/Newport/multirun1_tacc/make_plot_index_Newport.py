"""
Code to make plot indices for  multiple GeoClaw runs,
with a different dtopo file for each.
"""

import os,sys,shutil,glob
import numpy

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set


location = 'Newport'
show_gaugenos = [1007, 1022, 1041]

gaugemap_name = '%sGauges' % location
gaugemap_kml = '%s.kml' % gaugemap_name
gaugemap_jpg = '%s.jpg' % gaugemap_name
gaugemap_kml_path = root_dir + '/gauges/oregon_gauges_2024/%s' % gaugemap_kml
gaugemap_jpg_path = root_dir + '/gauges/oregon_gauges_2024/%s' % gaugemap_jpg

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
#geoclaw_plots = os.path.join(runs_dir, 'geoclaw_plots')


geoclaw_plots = './geoclaw_plots'  # for local version
#geoclaw_plots = './NewportPlots4hr'  # for local version

print('geoclaw_plots = ',geoclaw_plots)

print('Copying gaugemap files from\n    %s\n    %s' \
        % (gaugemap_kml_path, gaugemap_jpg_path))
shutil.copy(gaugemap_kml_path, geoclaw_plots)
shutil.copy(gaugemap_jpg_path, geoclaw_plots)

#cmd = 'cp ../%sGauges.* %s' % (location,geoclaw_plots)
#print(cmd)
#os.system(cmd)



# Select set of events to show plots for:


#CoPes Hub ground motions:

all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']

if 1:
    #models = all_models
    models = all_models[:2]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

events.sort()
events = events[:12]

instant = False
if instant:
    events = [e+'_instant' for e in events]
    rupture_type = 'Instant'




plotdirs = ['%s/_plots_%s' % (geoclaw_plots, event) for event in events]



def make_html_index_gauges(plotdir,event,gauges_plotdir):

    gauge_location_text = \
    """<h2> Gauge locations:</h2>
    Download <a href="../../%sGauges.kml">%sGauges.kml</a> and open in
    Google Earth.
    <p>
    <img src=../../%sGauges.jpg width=50%%>
    <p>
    &nbsp;
    <p>
    """ % (location,location,location)

    html_fname = os.path.join(plotdir,gauges_plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write(gauge_location_text)
        gauge_pngs = glob.glob(gauges_plotdir + '/*_Gauge*.png')
        gauge_pngs = glob.glob(os.path.join(plotdir, gauges_plotdir,
                                            '*_Gauge*.png'))
        gauge_pngs.sort()
        f.write('\n<p><h2>Gauge plots:</h2>\n')

        for k,gauge_png in enumerate(gauge_pngs):
            gauge_png = os.path.split(gauge_png)[-1]
            gaugename = os.path.splitext(gauge_png)[0]
            gaugeno = int(gaugename[-5:])
            f.write('<a href="%s">%i</a> &nbsp;&nbsp;\n' % (gauge_png, gaugeno))
            if numpy.mod(k+1,10) == 0:
                f.write('<br>\n')
        f.write('<p>\n')
    print('Created ',html_fname)


def make_html_index(plotdir,event,gauges_plotdir,fgmax_plotdir):
    """
    Make index for a single event.
    This function is called in loop below.
    """
    html_fname = os.path.join(plotdir,'index.html')
    with open(html_fname, 'w') as f:
        f.write('<html>\n<h1>%s</h1>\n' % event)
        f.write('\n<ul>\n<li><a href="%s">fgmax plots</a>\n' % fgmax_plotdir)
        f.write('<li><a href="%s/index.html">gauge plots</a>\n</ul>\n' % gauges_plotdir)
    print('Created ',html_fname)

    # make index for all gauges for this event:
    make_html_index_gauges(plotdir,event,gauges_plotdir)


# top_index show select plots from each event with links to plot index for each
top_index_fname = os.path.join(geoclaw_plots,'index.html')

gauge_location_text = \
"""<h2> Gauge locations:</h2>
Download <a href="%sGauges.kml">%sGauges.kml</a> and open in
Google Earth.
<p>
<img src=%sGauges.jpg width=50%%>
<p>
&nbsp;
<p>
""" % (location,location,location)

with open(top_index_fname, 'w') as top_index:

    #top_index.write('<html>\n<h1>Plots for %s</h1>\n<ul>\n' % location)
    top_index.write('<html>\n<h1>Inundation plots for %s</h1>\n' % location)
    top_index.write(
"""
Work in progress, supported in part by the NSF-funded
<a href="https://cascadiacopeshub.org/">Cascadia CoPes Hub</a>.
<p>
Tsunami inundation computed from new CoPes Hub CSZ ground motions
using the GeoClaw tsunami model.  The maximum flow depth and flow speed are
shown for a 2-hour simulation, both over the Yaquina Bay region and on 3
transects through this region.  Additional fgmax plots for this region are
linked as fgmax1 below each event,
while fgmax2 refers to the region to the east of longitude -124.
<p>
Synthetic gauge locations are shown below.  For each event, sample gauge
plots are shown at three gauges: <br>
<ul>
<li>1007 in the channel leading into Yaquina Bay, <br>
<li>1022 on shore near the OSU Marine Science Building, and <br>
<li>1041 where Yaquina River meets the Bay.<br>
</ul>
Additonal gauge plots and a report on maximum values seen at each gauge are
linked below each event.


""")

    top_index.write(gauge_location_text)

    for k in range(len(plotdirs)):
        plotdir = plotdirs[k]
        event = events[k]
        run_name = '%s_%s' % (location,event)

        gauges_plotdir = 'gauges'
        fgno = 1  # Yaquina Bay
        fgmax_plotdir = 'fgmax%s' % fgno

        make_html_index(plotdir,event,gauges_plotdir,fgmax_plotdir)

        # relative paths:
        plotdir = '_plots_%s' % event

        gauges_plotdir = '%s/%s' % (plotdir, gauges_plotdir)
        fgmax_plotdir = '%s/%s' % (plotdir, fgmax_plotdir)

        top_index.write('<h2>%s</h2>\n ' % event)

        # images to appear on index page:

        # fgmax plots
        top_index.write('<img src=%s/%s_fgmax%s_h_onshore.png width=28%%>\n' \
                % (fgmax_plotdir,run_name,fgno))
        top_index.write('<img src=%s/%s_fgmax%s_transects.png width=38%%>\n' \
                % (fgmax_plotdir,run_name,fgno))
        top_index.write('<img src=%s/%s_fgmax%s_speed.png width=28%%>\n<p>\n' \
                % (fgmax_plotdir,run_name,fgno))

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

            for fgno in [1,2]:
                fgmax_plotdir = fgmax_plotdir[:-1] + str(fgno)
                top_index.write('<li>  fgmax%s: <a href="%s/%s_fgmax%s_h_onshore.png">max h onshore</a>\n' \
                        % (fgno,fgmax_plotdir,run_name,fgno))
                top_index.write('&nbsp;&nbsp; <a href="%s/%s_fgmax%s_speed.png">max speed</a>\n' \
                        % (fgmax_plotdir,run_name,fgno))
                top_index.write('&nbsp;&nbsp; <a href="%s/%s_fgmax%s_transects.png">transects</a>\n' \
                        % (fgmax_plotdir,run_name,fgno))
                top_index.write('&nbsp;&nbsp; <a href="%s/%s_fgmax%s.kmz">kmz file</a>\n' \
                        % (fgmax_plotdir,run_name,fgno))

            top_index.write('<li> Gauges: <a href="%s/index.html">gauges</a>\n' \
                    % gauges_plotdir)
            top_index.write('&nbsp;&nbsp; <a href="%s/gauges_report_timesarrival.txt">gauges_report_timesarrival.txt</a>\n' \
                    % gauges_plotdir)
            top_index.write('&nbsp;&nbsp; <a href="%s/gauges_report_timesarrival.csv">gauges_report_timesarrival.csv</a>\n' \
                    % gauges_plotdir)



        top_index.write('</ul>\n')

print('Created', top_index_fname)
