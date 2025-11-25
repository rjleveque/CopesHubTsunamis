
from pylab import *
import numpy
import clawpack.pyclaw.gauges as gauges
import clawpack.amrclaw.data as amrclaw

from scipy.interpolate import interp1d
import os
import netCDF4
import xarray
import time as time_module


def make_all_gauges_nc(location, events, outdirs, gaugenos,
                       nc_fname=None, dt=5):

    drytol = 1e-3  # set u=v=0 where h < drytol when computing from hu,hv

    # quantities of interest:
    gauge_qois = ['h','u','v','eta','level']
    iqois = range(len(gauge_qois))  # integer indices 0,1,2,3
    qois = gauge_qois

    # netcdf data type:
    datatype = 'f4'  # 'f4' for 4-byte float32, 'f8' for 8-byte float64

    # times:
    tmin = inf
    tmax = -inf

    missing = []
    for k,event in enumerate(events):
        outdir = f'{outdirs}/_output_{event}'
        if not os.path.isdir(outdir):
            print(f'*** Skipping event {event},')
            print('***   no such outdir = ', outdir)
            missing.append(event)
        else:
            for j,gaugeno in enumerate(gaugenos):
                gauge = gauges.GaugeSolution(gaugeno, outdir)
                t = gauge.t
                tmin = min(tmin, t.min())
                tmax = max(tmax, t.max())

    if len(missing) == len(events):
        print('*** all events missing')
        return


    print('minimum time found at any gauge = %.2f seconds' % tmin)
    print('maximum time found at any gauge = %.2f seconds' % tmax)
    tg = arange(tmin, tmax+.1*dt, dt)
    print(f'Gauge time series will use tg of length {len(tg)} with dt = {dt}')

    if nc_fname is None:
        print('Must specify nc_fname to create netCDF file')
        return

    gauge_results = nan*zeros((len(tg),len(gaugenos),len(qois),len(events)),
                          dtype=float)

    xg = nan*zeros(len(gaugenos))
    yg = nan*zeros(len(gaugenos))


    for k,event in enumerate(events):
        outdir = '%s/_output_%s'  % (outdirs,event)
        for j,gaugeno in enumerate(gaugenos):
            gauge = gauges.GaugeSolution(gaugeno, outdir)
            t = gauge.t

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

            level = gauge.level
            gaugefcn = interp1d(t, level, kind='linear', bounds_error=False)
            gauge_results[:,j,4,k] = gaugefcn(tg)


    with netCDF4.Dataset(nc_fname, 'w') as rootgrp:

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

    print(f'Created {nc_fname}')
