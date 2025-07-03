"""
Tools under development for the Cascadia CoPes Hub project
(https://cascadiacopeshub.org/) supported by NSF.

"""

from pylab import *
import xarray
from scipy.interpolate import RegularGridInterpolator


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
    """

    @property
    def X(self):
        r"""Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates(mask=False)
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
            self.generate_2d_coordinates(mask=False)
        return self._Y
    @Y.setter
    def Y(self, value):
        self._extent = None
        self._Y = value
        self._y = nan

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

    The netCDF file can be generated using
        user_tools/make_allgauges_nc.py

    Time series for some quantities of interest (qoi) have been saved,
    typically ['h' 'u' 'v' 'eta'], where h is the water depth and
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
            4-dimensional array with all time series values


    """

    import xarray

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
