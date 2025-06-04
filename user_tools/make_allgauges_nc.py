
from pylab import *
import numpy
import clawpack.pyclaw.gauges as gauges
import clawpack.amrclaw.data as amrclaw

from scipy.interpolate import interp1d
import os
import netCDF4
import xarray
import time as time_module


dry_run = True

drytol = 1e-3  # set u=v=0 where h < drytol when computing from hu,hv

datatype = 'f4'  # 'f4' for 4-byte float32, 'f8' for 8-byte float64


outdirs = '/gscratch/tsunami/rjl/CopesHubTsunamis/geoclaw_runs/sites/seaside/multirun2/geoclaw_outputs'

location = 'Seaside'

nc_file = 'allgauges_%s.nc' % location
print('nc_file = ',nc_file)

# make a netcdf file for this one event with these gauges:
gaugenos = range(1001,1068,1)


# specify events...

all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


models = all_models
events = ['%s-deep' % model for model in models] \
       + ['%s-middle' % model for model in models] \
       + ['%s-shallow' % model for model in models]

events.sort()

#events = events[9:]  # only "random" events
events = [e for e in events if 'random-str10' in e]  # only "random-str10"

# or test on a few:
#events = ['buried-random-mur13-deep']

# quantities of interest:
gauge_qois = ['h','u','v','eta']
iqois = range(len(gauge_qois))  # integer indices 0,1,2,3
qois = gauge_qois

# times:
t0 = 0.
tfinal = 5400.
dt = 5.
tg = arange(t0, tfinal+.1*dt, dt)

gauge_results = empty((len(tg),len(gaugenos),len(qois),len(events)),
                      dtype=float)

xg = nan*zeros(len(gaugenos))
yg = nan*zeros(len(gaugenos))

tmin = inf
tmax = -inf

if dry_run:
    print('dry_run: will check that output directories exist for these events:')
    print(events)
missing_directories = False

for k,event in enumerate(events):
    outdir = '%s/_output_%s'  % (outdirs,event)
    if os.path.isdir(outdir):
        print('Found directory: ',outdir)
    else:
        print('Missing directory: ',outdir)
        missing_directories = True

if missing_directories:
    raise ValueError('*** aborting due to missing directories')

for k,event in enumerate(events):
    outdir = '%s/_output_%s'  % (outdirs,event)
    for j,gaugeno in enumerate(gaugenos):
        gauge = gauges.GaugeSolution(gaugeno, outdir)
        t = gauge.t
        tmin = min(tmin, t.min())
        tmax = max(tmax, t.max())

        x,y = gauge.location
        if isnan(xg[j]):
            # first event, set x,y for this gauge:
            xg[j] = x
            yg[j] = y
        if x != xg[j] or y != yg[j]:
            print('*** x or y do not match for gaugeno %i, event %s' \
                    % (gaugeno,event))
            print('*** x,y = %g,%g, expected %g,%g' \
                    % (x,y,xg[j],yg[j]))

        if dry_run:
            break

        # set the data:
        h = gauge.q[0,:]
        gaugefcn = interp1d(t, h, kind='linear', bounds_error=False)
        gauge_results[:,j,0,k] = gaugefcn(tg)

        u = numpy.divide(gauge.q[1,:], h,
                          out=numpy.zeros(h.shape, dtype=float), \
                          where=(h>drytol))
        gaugefcn = interp1d(t, u, kind='linear', bounds_error=False)
        gauge_results[:,j,1,k] = gaugefcn(tg)

        v = numpy.divide(gauge.q[2,:], h,
                          out=numpy.zeros(h.shape, dtype=float), \
                          where=(h>drytol))
        gaugefcn = interp1d(t, v, kind='linear', bounds_error=False)
        gauge_results[:,j,2,k] = gaugefcn(tg)

        eta = gauge.q[-1,:]
        gaugefcn = interp1d(t, eta, kind='linear', bounds_error=False)
        gauge_results[:,j,3,k] = gaugefcn(tg)

    if dry_run:
        break

print('minimum time found at any gauge = %.2f seconds' % tmin)
print('maximum time found at any gauge = %.2f seconds' % tmax)

if dry_run:
    print('Rerun with dry_run==False to make netCDF file')

if not dry_run:
    with netCDF4.Dataset(nc_file, 'w') as rootgrp:

        rootgrp.description = 'All gauge output for %s'  % location

        rootgrp.history = 'Created ' + time_module.ctime(time_module.time())
        rootgrp.outdir = os.path.abspath(outdir)
        #rootgrp.run_finished = run_finished

        #gauge = rootgrp.createGroup(gauge_name)
        time = rootgrp.createDimension('time', len(tg))
        gaugeno = rootgrp.createDimension('gaugeno', len(gaugenos))
        #iqoi = rootgrp.createDimension('iqoi', len(iqois))
        qoi = rootgrp.createDimension('qoi', len(qois))
        event = rootgrp.createDimension('event', len(events))

        gaugeno = rootgrp.createVariable('gaugeno','i4',('gaugeno',))
        gaugeno[:] = gaugenos

        event = rootgrp.createVariable('event','str',('event',))
        for j,e in enumerate(events):
            event[j] = e

        qoi = rootgrp.createVariable('qoi','str',('qoi',))
        for j,q in enumerate(qois):
            qoi[j] = q

        x = rootgrp.createVariable('x','f8',('gaugeno',))
        x[:] = xg
        x.units = "longitude_east"

        y = rootgrp.createVariable('y','f8',('gaugeno',))
        y[:] = yg
        y.units = "latitude_north"

        times = rootgrp.createVariable('time',datatype,('time',))
        times[:] = tg
        times.units = "seconds"

        gauge_vals = rootgrp.createVariable('gauge_vals',datatype,
                            ('time','gaugeno','qoi','event'))
        gauge_vals[:,:,:,:] = gauge_results

    print('Created %s' % nc_file)
