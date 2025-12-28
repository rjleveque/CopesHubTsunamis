import os
import numpy as np

# top level directory for this project:
root_dir     = os.environ['CHT']

rundir = os.getcwd()

if '/home/ptha' in rundir:
    # on cubano:
    topo_dir = '/home/ptha/topo/topofiles'
    dtopo_dir = '/home/ptha/dtopo/dtopofiles'
elif 'mmfs1' in rundir:
    # on hyak:
    topo_dir = os.path.abspath('/gscratch/tsunami/topo/topofiles/')
    dtopo_dir = os.path.abspath('/gscratch/tsunami/dtopo/dtopofiles/')
else:
    # on laptop  # ADJUST!
    topo_dir =  root_dir + '/topo/topofiles'
    dtopo_dir = root_dir + '/dtopo/CSZ_groundmotions/dtopofiles'

# for dtopofiles in this project root_dir use this instead of above:
# except for CSZ dtopo which is in dtopo_dir above
# dtopo_dir = root_dir + '/dtopo/dtopofiles'

# Set out_dir and plot_dir, need for making plots.
# Should be consistent with OUTDIR and PLOTDIR specified in Makefile.

if 'mmfs1' in rundir:
    # on hyak output and plots on scratch disk:
    scrdir = rundir.replace('mmfs1/home','gscratch/tsunami')
    out_dir = os.path.join(scrdir, '_output')
else:
    # standard location:
    out_dir = os.path.join(rundir, '_output')
plot_dir = out_dir.replace('_output','_plots')


loc = 'CHT'
event = 'buried-locking-str10-deep_instant'
dtopofile = dtopo_dir + '/buried-locking-str10-deep_instant.dtt3'

tfinal = 1*3600.
#tfinal = 2700.
#Doing a 15 minute test run just to get the subsidence at gauges
#and initial B0, etc.
#tfinal = 15*60.
#dtout = 10*60.
#output_times = [0.,60.] + list(np.arange(dtout, tfinal+1, dtout))
output_times = [0.,60.,tfinal]

##Need to fix this, but program doesnt use it
fgmax_extent=[-124.1825,-124.135,46.975,46.9875]

