from pylab import *
import os,sys,glob
import datetime

plotdir = '_plots'
fname_html = '%s/index.html' % plotdir

with open(fname_html,'w') as f:
    f.write('<html>\n<h1>Gauge plots for event ASCE_SIFT_Region2</h1>\n')
    timenow = datetime.datetime.today().strftime('%Y-%m-%d at %H:%M:%S')
    f.write('<h3>Generated %s</h3>\n' % timenow)
    
    f.write("""<h2>Tsunami propagation</h2>\n
This <a href="geoclaw_ASCE_SIFT_Region2.mp4">animation</a>
shows the tsunami propagation.
<p>
This page shows a study of how the resolution of the topography DEM and/or
the computational grids used affects the maximum wave amplitude (elevation
above MSL) when the same event is used.  Later studies will compare
different events or time-dependent vs. instantaneous ruptures.\n""")


    f.write('<h2>Maximum amplitude at each gauge for same event:</h2>\n')
    f.write('<p><img src="rotated_map_gauge_comparisons.png" width=1200>\n')
    
    f.write('<h2>Zoomed view around OR and WA:</h2>\n')
    f.write('<p><img src="rotated_map_gauge_comparisons_zoom.png" width=1200>\n')

    f.write('<h2>Time series for select ASCE gauges:</h2>\n')
    
    for gaugeno in [161,125,182]:
        f.write('<p><img src="gauge%s_comparison.png" width=1200>\n' \
                    % str(gaugeno).zfill(5))

    f.write('<h2>Time series for other ASCE gauges:</h2>\n')
    
    f.write('Gauge number ...\n')
    for gaugeno in range(1,651,1):
        f.write('<a href="gauge%s_comparison.png">%s</a>' \
                    % (str(gaugeno).zfill(5), gaugeno) + '&nbsp ... &nbsp')                  
    f.write('\n\n</html>\n')
