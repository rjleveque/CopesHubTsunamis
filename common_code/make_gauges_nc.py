
"""
Sample code to make a single .nc file with gauge outputs for all gauges
at one location, from all events in a specified list.

Run this code from the multirun directory
(from which runclaw_makeplots_dtopos.py was run)

Before running, adjust to set:
    events to the desired set of events
    location to the string that will appear in the gauges .nc file
"""
import os,sys,glob
import pathlib
from clawpack.visclaw import gaugetools
from clawpack.clawutil.util import fullpath_import

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

try:
    CHTuser = os.environ['CHTuser']
except:
    raise Exception("*** Must first set CHTuser environment variable")

CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')

location = 'Aberdeen' # to appear in .nc file name



# location for big files for different computer environments:
this_dir = os.getcwd()
HOME = os.environ['HOME']

if '/home1' in this_dir:
    computer = 'tacc'
    #scratch_dir = this_dir.replace('/home1', '/scratch')
    try:
        SCRATCH = os.environ['SCRATCH']
        scratch_dir = this_dir.replace(HOME, SCRATCH)
    except:
        scratch_dir = this_dir  # if $SCRATCH not set

# geoclaw_outputs should be full path to the directory that contains
# subdirectories such as _output_BL10D for individual event runs:
geoclaw_outputs = f'{scratch_dir}/geoclaw_outputs'


if 0:
    # specify which events to include:

    events = CHTtools.all_events()
    events = events[:4]

    #events = ['BL10D', 'BL13D', 'BL13M']

    instant = True
    if instant:
        events = [f'{e}_instant' for e in events]

else:
    # include all events found in geoclaw_outputs:
    outdirs = glob.glob(f'{geoclaw_outputs}/_output_*')
    events = []
    for outdir in outdirs:
        p = pathlib.Path(outdir)
        name = p.stem
        event = name.replace('_output_','')
        events.append(event)

print(f'Using {len(events)} events: ',events)
print(f'      from {geoclaw_outputs}')


# assume all events have same set of gauges, read from first outdir:
outdir0 = f'{geoclaw_outputs}/_output_{events[0]}'
setgauges = gaugetools.read_setgauges(outdir0)
gaugenos = setgauges.gauge_numbers
print('Found gaugenos = ',gaugenos)

# or set list of gauges to include explicitly:
#gaugenos = [101,102]

nc_fname = f'{location}_gauges_4events.nc'

CHTtools.make_all_gauges_nc(location, events, geoclaw_outputs, gaugenos,
                            nc_fname, dt=5)
