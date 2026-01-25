
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

events = CHTtools.all_events()
events = events[:4]

#events = ['BL10D', 'BL13D', 'BL13M']

instant = True
if instant:
    events = [f'{e}_instant' for e in events]


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

location = 'Aberdeen' # to appear in file name

setgauges = gaugetools.read_setgauges(outdirs[0])
gaugenos = setgauges.gauge_numbers
print('Found gaugenos = ',gaugenos)

nc_fname = f'{location}_gauges_test.nc'
#nc_fname = None

CHTtools.make_all_gauges_nc(location, events, outdirs, gaugenos,
                                  nc_fname, dt=5)
