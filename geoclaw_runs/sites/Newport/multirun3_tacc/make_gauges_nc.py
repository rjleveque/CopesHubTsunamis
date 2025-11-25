
import os,sys,glob

from clawpack.visclaw import gaugetools
from clawpack.clawutil.util import fullpath_import

try:
    CHT = os.environ['CHT']
except:
    raise Exception("*** Must first set CHT environment variable")

user_tools_dir = os.path.join(CHT, 'user_tools')
CHTtools = fullpath_import(f'{user_tools_dir}/CHTtools.py')
CHT_gaugetools = fullpath_import(f'{user_tools_dir}/CHT_gaugetools.py')

#events = CHTtools.all_events()

events = ['BL10D', 'BL13D', 'BL13M']

outdirs = './geoclaw_outputs'

location = 'Newport'

setgauges = gaugetools.read_setgauges(outdirs[0])
gaugenos = setgauges.gauge_numbers
print('Found gaugenos = ',gaugenos)

nc_fname = f'{location}_gauges_test.nc'
#nc_fname = None

CHT_gaugetools.make_all_gauges_nc(location, events, outdirs, gaugenos,
                                  nc_fname, dt=5)
