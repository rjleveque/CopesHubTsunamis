"""
Tools under development for the Cascadia CoPes Hub project
(https://cascadiacopeshub.org/) supported by NSF.

"""

from pylab import *
import xarray
from scipy.interpolate import RegularGridInterpolator

def all_events():
    depths = ['D','M','S']

    # buried_locking events:
    all_events = [f'BL10{depth}' for depth in depths] \
               + [f'BL13{depth}' for depth in depths] \
               + [f'BL16{depth}' for depth in depths] \

    # add random events:
    all_events += [e.replace('L','R') for e in all_events]

    # add ft events:
    all_events += [e.replace('B','F') for e in all_events]

    return all_events


def shortname(event):
    """
    convert original event name like ft-locking-mur13-deep
    to new short name like FL13D
    """

    if len(event) == 5:
        # assume it is already a short name:
        return event

    tokens = event.replace('_','-').split('-')
    newname = ''
    if tokens[0] == 'buried':
        newname += 'B'
    elif tokens[0] == 'ft':
        newname += 'F'
    if tokens[1] == 'locking':
        newname += 'L'
    elif tokens[1] == 'random':
        newname += 'R'
    newname += tokens[2][-2:]  # year
    if tokens[3] == 'deep':
        newname += 'D'
    elif tokens[3] == 'middle':
        newname += 'M'
    elif tokens[3] == 'shallow':
        newname += 'S'

    if len(newname) != 5:
        print(f'*** problem converting {event}, newname = {newname}')

    return newname

def name_conversions():
    models = \
       ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
        'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

    events = events + [e.replace('buried','ft') for e in events]
    events.sort()
    for event in events:
        #print(f'  {event.ljust(30)} --> {shortname(event)}')
        print(f'  {shortname(event)} = {event}')


class TsunamiModelResults(object):

    """
    General class for encapsulating tsunami modeling results from either
    GeoClaw or MOST, using the same names for qoi's (quantities of interest)
    such as depth, eta (surface), B (topo/bathy), etc. and the same
    index ordering corresponding to (x,y) for 2D arrays like B and
    (x,y,t) for 3D arrays like depth.

    Functions are included for reading in results on a fixed spatial grid
    at a sequence of times from a netCDF file produced from the tsunami
    simulation, and for creating interpolating functions that allow the
    evaluation of a qoi at arbitrary points in space and time.

    Still under development and may change.

    Methods defined for this class:

    generate_2d_coordinates:
        After 1D arrays x,y have been defined as attributes, e.g. by
        load_results(), this is used to create 2D versions X,Y.
        Invoked automatically when X,Y attributes are used (by @property).

    load_results:
        load tsunami simulation from a netCDF file using either 'geoclaw'
        or 'most' format. Sets attributes such as `depth` as 3D arrays
        with indices [i,j,k] refering to (x[i],y[j],t[k]) and
        `B0` (inital topography) as 2D array with indices [i,j].

    make_interp_fcn:
        convert an attribute such as `depth` into a function that can be
        evaluated at any (x,y,t) within the domain and time range in the
        dataset. (Returns nan at points outside).

    """

    @property
    def X(self):
        r"""Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates()
        return self._X
    @X.setter
    def X(self, value):
        self._extent = None
        self._X = value
        self._x = nan

    @property
    def Y(self):
        r"""Two dimensional coordinate array in y direction."""
        if self._Y is None:
            self.generate_2d_coordinates()
        return self._Y
    @Y.setter
    def Y(self, value):
        self._extent = None
        self._Y = value
        self._y = nan

    @property
    def extent(self):
        r"""Extent of domain"""
        if self._extent is None:
            self.generate_2d_coordinates()
        return self._extent


    def __init__(self, x=None, y=None, t=None, ncfile=None, format=None):
        self.x = x
        self.y = y
        self.t = t
        self.ncfile = ncfile
        self.format = format

        # possible qois (quantities of interest):
        self.depth = None
        self.u = None
        self.v = None
        self.speed = None
        self.momFlux = None
        self.B = None
        self.B0 = None
        self.Bfinal = None
        self.qois = []  # will be set to qois found in ncfile

        # 2D arrays only calculated when needed (using @property):
        self._X = None
        self._Y = None
        self._extent = None


    def generate_2d_coordinates(self):
        self._X, self._Y = meshgrid(self.x, self.y, indexing='ij')
        self._extent = [self.x.min(), self.x.max(),
                        self.y.min(), self.y.max()]


    def load_results(self, qois='all'):

        with xarray.open_dataset(self.ncfile, decode_timedelta=False) as ncdata:

            if self.format == 'most':
                xdim = 'x'
                ydim = 'y'
                mapping = {'depth':'depth', 'eta':'h', 'speed':'speed',
                           'u':'u', 'v':'v', 'B0':'bathy', 'Bfinal':'Bfinal',
                           'B':'B', 'momFlux':'momFlux'}

            elif self.format == 'geoclaw':
                xdim = 'lon'
                ydim = 'lat'
                mapping = {'depth':'h', 'eta':'eta', 'speed':'s',
                           'u':'u', 'v':'v', 'B0':'B0', 'Bfinal':'Bfinal',
                           'B':'B', 'momFlux':'momflux'}
            else:
                raise InputError('Unrecognized format: ',format)

            self.x = asarray(ncdata.variables['lon'])
            self.y = asarray(ncdata.variables['lat'])
            self.t = asarray(ncdata.variables['time'])

            if qois == 'all':
                qois = ['depth','eta','u','v','B']

            qois_found = []

            for qoi in qois:
                try:
                    vals = ncdata.variables[mapping[qoi]]

                    # reorganize multidimensional arrays if necessary, so that
                    #    [i,j,k] index corresponds to (x,y,t)
                    # and then convert to numpy arrays:
                    vals = asarray(vals.transpose(xdim,ydim,'time'))
                    setattr(self, qoi, vals)
                    print('qoi %s(x,y,t) set with shape = %s' \
                              % (qoi,vals.shape))
                    qois_found.append(qoi)
                except:
                    print('qoi %s not found in %s' % (qoi,self.ncfile))

            for qoi in ['B0', 'Bfinal']:
                try:
                    vals = ncdata.variables[mapping[qoi]]

                    # reorganize multidimensional arrays if necessary, so that
                    #    [i,j] index corresponds to (x,y)
                    # and then convert to numpy arrays:
                    vals = asarray(vals.transpose(xdim,ydim))
                    setattr(self, qoi, vals)
                    print('qoi %s(x,y) set with shape = %s' \
                              % (qoi,vals.shape))
                    qois_found.append(qoi)
                except:
                    print('qoi %s not found in %s' % (qoi,self.ncfile))

        self.qois = qois_found

    def make_interp_fcn(self, qoi):
        """
        make interpolator to evaluate at arbitrary (x,y) or (x,y,t)
        """

        method = 'linear'
        bounds_error = False
        fill_value = nan

        vals = getattr(self, qoi)

        if vals.ndim == 2:
            # e.g. B0
            qoi_fcn = RegularGridInterpolator((self.x,self.y), vals,
                                              method,bounds_error,fill_value)
        elif vals.ndim == 3:
            qoi_fcn = RegularGridInterpolator((self.x,self.y,self.t), vals,
                                              method,bounds_error,fill_value)

        return qoi_fcn

    def make_all_qoi_fcns(self):
        for qoi in self.qois:
            qoi_fcn = self.make_interp_fcn(qoi)
            name_qoi_fcn = '%s_fcn' % qoi
            setattr(self, name_qoi_fcn, qoi_fcn)
            print('created function %s' % name_qoi_fcn)

    # end class definition


def load_tsunami_model(ncfile, format='geoclaw'):

    with xarray.open_dataset(ncfile, decode_timedelta=False) as ncdata:

        t = asarray(ncdata.variables['time'])

        if format == 'most':
            xdim = 'x'
            ydim = 'y'
            x = asarray(ncdata.variables['lon'])
            y = asarray(ncdata.variables['lat'])

        elif format == 'geoclaw':
            xdim = 'lon'
            ydim = 'lat'
            x = asarray(ncdata.variables['lon'])
            y = asarray(ncdata.variables['lat'])

        else:
            raise InputError('Unrecognized format: ',format)


    model = TsunamiModelResults(x,y,t,ncfile,format)
    model.load_results()

    return model



def read_allgauges_nc(ncfile):

    """
    Read all gauge data for a particular coastal site from a netCDF file
    that contains time series at all gauges pre-selected in this region,
    for some set of earthquake events.

    The netCDF file can be generated using the function
    `make_allgauges_nc` defined below.

    Time series for some quantities of interest (qoi) have been saved,
    typically ['h' 'u' 'v' 'eta' 'level'], where h is the water depth and
    eta is the surface elevation relative to the vertical datum of the
    topography file (e.g. MHW).

    :Input:
        - ncfile : str
            path to the netCDF file

    :Output:
        - gauge_x : xarray.DataArray
            longitude of each gauge, indexed by gaugeno
        - gauge_y : xarray.DataArray
            latitude of each gauge, indexed by gaugeno
        - gauge_t : numpy.ndarray
            times for the time series, in seconds
            for convenience in plotting, with values from
            gauge_vals.coords.time
        - gauge_vals : xarray.DataArray
            5-dimensional array with all time series values


    """

    import xarray

    print(f'Trying to read {ncfile}  ...')

    with xarray.open_dataset(ncfile, decode_timedelta=False) as ncdata:
        print(ncdata.description)
        gaugeno = ncdata.variables['gaugeno']
        x = ncdata.variables['x']
        y = ncdata.variables['y']
        time = ncdata.variables['time']
        event = ncdata.variables['event']
        qoi = ncdata.variables['qoi']
        gauge_vals = ncdata.variables['gauge_vals']

        # create xarrays to pass back:
        dims = ('time','gaugeno','qoi','event')
        coords = (time, gaugeno, qoi, event)
        gauge_vals = xarray.DataArray(gauge_vals, coords, dims)

        gauge_x = xarray.DataArray(x, (gaugeno,), ('gaugeno',))
        gauge_y = xarray.DataArray(y, (gaugeno,), ('gaugeno',))
        gauge_t = time.data  # simple numpy.ndarray for convenience

    print('Loaded %i times from t = %.1f to %.1f sec at %i gauges' \
            % (len(time),time[0],time[-1],len(x)))
    print('   for %i events, with qois %s' % (len(event), qoi.data))
    print('print gauge_vals.coords for more info')

    return gauge_x, gauge_y, gauge_t, gauge_vals


def make_all_gauges_nc(location, events, outdirs, gaugenos,
                       nc_fname=None, dt=5):
    """
    Make a netCDF file containing all specified `gaugnos` for all `events`.
    Assumes `outdirs` has subdirectories with names like `_output_BL13D`.
    Uses pw linear interpolation with time increment `dt` for output.
    """

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
