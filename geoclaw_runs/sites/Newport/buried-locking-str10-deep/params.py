
import os

location = 'Newport'  # for naming plots and mp4 file

rundir = os.getcwd()
print('rundir = %s' % rundir)

# for new ground motions:
event = os.path.split(os.getcwd())[-1]
instant = ('instant' in rundir)  # is this instantaneous uplift?

print('event = %s, instant = %s' % (event,instant))
