#!/usr/bin/env python
# coding: utf-8

# # Adapted from OkadaTest.ipynb 241023
#
# Read in a fault model from the coarsened triangulation provided by Jey and apply the Okada model from GeoClaw to obtain the static seafloor deformation at the final time.
#
# (Could be adapted to compute time-dependendent kinetic rupture based on applying Okada at each specified time to the subfaults that have ruptured up to that time, but that would require `rupture time` and `rise_time` for each subfault, which are not included currently in Jey's coarsened model.)
#


from pylab import *
#from clawpack.geoclaw import dtopotools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import copy
import os,sys
from multiprocessing import current_process


# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.join(root_dir, 'common_code'))
from multip_tools import run_many_cases_pool
import dtopotools  # from common_code, refactored use of subfault.dtopo

dry_run = False  # if True, just print event names
test_subset = False

mu = 30e9  # rigidity = shear modulus (in Pascals)
# Note: the value of mu doesn't matter in Okada solution, only for computing Mw
# Okada depends only on Poisson ratio which is set to 0.25 by default.


rupture_type = 'static'
    # 'static' ==> instant displacement at t=0
    # 'kinematic' ==> apply Okada to each subfault but only include dz in
    #                 dtopo.dZ for times after subfault.rupture_time
    #                 times will be every 10 seconds to end of rupture



# ## Read in the fault geometry:
#
# And specify a `dtopotools.Fault` object with this geometry.


#datadir = '/Users/rjl/D/JeysCode/Audreys_source_models_interpolation_09182024/'
#triangles = loadtxt(datadir+'David_cas_fine_mesh.tri')  # list of triangles, 3 vertices
#vertices = loadtxt(datadir+'David_cas_fine_mesh.ned')  # list of vertices, lon-lat of each

datadir = 'new_coarsened_source_random_mur13_deep/'
triangles = loadtxt(datadir+'JK_cas_fine_mesh_updip_modified.tri')  # list of triangles, 3 vertices
vertices = loadtxt(datadir+'JK_cas_fine_mesh_updip_modified.ned')  # list of vertices, lon-lat of each

vertices[:,1] = vertices[:,1] # - 360.  # shift to longitude W
vertices[:,3] = -1000*vertices[:,3]   # convert depth to positive depth in meters


# Create Fault with geometry, no slip yet:

fault0 = dtopotools.Fault(coordinate_specification='triangular')
fault0.rupture_type = rupture_type
fault0.subfaults = []


# Grid for dtopo vertical deformations:
dx = dy = 30/3600.  # spatial resolution for dtopo file
#x,y = fault0.create_dtopo_xy(dx=dx)  # choose automatically
# use same x,y as in dtopo files made from ground motions:
x = arange(-128.5, -122.4999, dx)
y = arange(40,50.0001,dy)

print('Will create dtopo on arrays of shape %i by %i ...' % (len(x),len(y)))


nsubfaults = triangles.shape[0]

for j in range(nsubfaults):
    subfault0 = dtopotools.SubFault()
    jv = int(triangles[j,1]) - 1
    node1 = vertices[jv,1:4]
    jv = int(triangles[j,2]) - 1
    node2 = vertices[jv,1:4]
    jv = int(triangles[j,3]) - 1
    node3 = vertices[jv,1:4]
    node_list = [node1, node2, node3]
    subfault0.set_corners(node_list,projection_zone='10')
    fault0.subfaults.append(subfault0)

print('Set up fault0 model with %i subfaults, without yet specifying slip for particular event' \
            % nsubfaults)


# ## Check that orientation of triangles are all correct:

numpos = 0.
for s in fault0.subfaults:
    c = array(s.corners)[:3,:2]
    A = vstack([c[:,0], c[:,1], array([1,1,1])]).T
    detA = np.linalg.det(A)
    if detA > 0:
        numpos += 1
if numpos > 0:
    print('*** Warning, %i of the %i subfaults have counterclockwise orientation' \
            % (numpos, nsubfaults))
else:
    print('All subfault triangles have clockwise orientation')


if 0:
    # ## Plot triangulation:

    fig = plt.figure(figsize=(15,10))
    #ax = fig.add_subplot(121, projection='3d')
    ax = fig.add_axes([.05,.05,.9,.9], projection='3d')
    for s in fault0.subfaults:
        c = s.corners
        c.append(c[0])
        c = np.array(c)
        ax.plot(c[:,0],c[:,1],-c[:,2]/1000.,color='b')
    ax.view_init(10,60)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth (km)')
    ax.set_title('Triangular subfaults')

    #ax = fig.add_subplot(122)
    #ax = fig.add_axes([.75,.05,.2,.9])
    fig = figure(figsize=(10,15))
    ax = axes()
    for s in fault0.subfaults:
        c = s.corners
        c.append(c[0])
        c = np.array(c)
        ax.plot(c[:,0],c[:,1], 'b')
    ax.set_aspect(1./np.cos(45*np.pi/180.))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Plan view');


# ## Functions for setting up an event and applying Okada:


def set_slip(fault0, event):

    # start with deepcopy of fault0 to make copies of each subfault
    # in fault0.subfaults, so that each subfault.slip can be modified
    # differently for each event:
    fault = copy.deepcopy(fault0)

    if 0:
        # alternative to deep copy:
        fault = dtopotools.Fault(coordinate_specification='triangular')
        fault.subfaults = []
        for s in fault0.subfaults:
            fault.subfaults.append(copy.copy(s))  # must be copy, not orig!

    #event_dir = datadir + 'output_files_combined'  # now includes kinematics
    event_dir = datadir

    # switch to Jey's notation:
    event1 = event.replace('buried-','')
    event_jey = event1.replace('-','_')
    print('+++ in set_slip, event = %s' % event)
    print('+++ in set_slip, event_jey = %s' % event_jey)

    if 0:
        # don't need slip in dip and strike directions separately, only magnitude read below
        fname = 'dip_slip_resampled_source_saved_%s.out' % event_jey
        dip_slip = loadtxt(os.path.join(event_dir, fname))
        fname = 'strike_slip_resampled_source_saved_%s.out' % event_jey
        strike_slip = loadtxt(os.path.join(event_dir, fname))

    fname = 'mag_slip_resampled_source_saved_%s.out' % event_jey
    mag_slip = loadtxt(os.path.join(event_dir, fname))[:,2]
    mag_slip = where(isnan(mag_slip), 0, mag_slip)  # replace nan values

    fname = 'rake_resampled_source_saved_%s.out' % event_jey
    rake = loadtxt(os.path.join(event_dir, fname))[:,2]
    rake = where(isnan(rake), 0, rake)  # replace nan values

    if rupture_type == 'kinematic':
        fname = 'rupt_time_resampled_source_saved_%s.out' % event_jey
        rupt_time = loadtxt(os.path.join(event_dir, fname))[:,2]
        rupt_time = where(isnan(rupt_time), 0, rupt_time)  # replace nan values
        fname = 'rise_time_resampled_source_saved_%s.out' % event_jey
        rise_time = loadtxt(os.path.join(event_dir, fname))[:,2]
        rise_time = where(isnan(rise_time), 0, rise_time)  # replace nan values
        rupt_tfinal = (rupt_time + rise_time).max()
        fault.dtopo_times = arange(0, rupt_tfinal+10, 10)  # every 10 seconds
    else:
        fault.dtopo_times = [0.]  # only compute static displacement (instant)

    print('+++ fault.rupture_type = %s, fault.dtopo_times = %s' \
                % (fault.rupture_type, fault.dtopo_times))

    for j in range(nsubfaults):
        subfault = fault.subfaults[j]
        subfault.rake = rake[j]
        subfault.slip = mag_slip[j]
        subfault.mu = mu
        if rupture_type == 'kinematic':
            subfault.rupture_time = rupt_time[j]
            subfault.rise_time = rise_time[j]
            subfault.rise_shape = 'quadratic'  # correct?
            subfault.rise_time_starting = None # correct?

    slips = array([s.slip for s in fault.subfaults])
    print('+++ max slip = %.2fm' % slips.max())
    print('In set_slip, created fault with Mw = %.2f' % fault.Mw())

    # save new name of event in fault object

    if rupture_type == 'static':
        fault.event = event + '_okada_instant'
    else:
        fault.event = event + '_okada_kinematic'

    print('New event name: %s' % fault.event)
    print('+++ fault.dtopo_times = ',fault.dtopo_times)
    return fault


def make_dtopo(fault):

    # Apply Okada to all subfaults to create deformation dtopo.dZ:
    dtopo = fault.create_dtopography(x,y,times=fault.dtopo_times,verbose=100)

    fname = '%s.dtt3' % fault.event
    dtopo.write(fname, dtopo_type=3)
    print('Created %s with %s displacement at %i times' \
            % (fname, fault.rupture_type, len(dtopo.times)))
    return dtopo


def plot_slip_final_dtopo(fault, dtopo):

    """
    Plot dtopo (at final time for kinematic rupture).
    """

    fig,(ax0,ax1) = plt.subplots(ncols=2,nrows=1,figsize=(12,6))
    fault.plot_subfaults(axes=ax0,slip_color=True,plot_box=False);
    ax0.set_title('Slip on Fault');

    X = dtopo.X; Y = dtopo.Y; dZ_at_t = dtopo.dZ_at_t
    dz_max = dtopo.dZ.max()
    tfinal = dtopo.times[-1] + 1  # 1 second after final dZ
    dtopotools.plot_dZ_colors(X,Y,dZ_at_t(tfinal),axes=ax1,
                              cmax_dZ = dz_max,
                              dZ_interval=200, add_colorbar=True);
    ax1.set_title('Seafloor deformation (static Okada)');

    fname = '%s.png' % fault.event
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)


def plot_slip_vs_seismic_instant(fault, dtopo):

    # ## Compare to final displacement from seismic simulation:

    dtopodir = '../dtopofiles/'
    dtopo_instant_fname = dtopodir + '%s.dtt3' % fault.event.replace('_okada','')
    print('Comparing to %s' % dtopo_instant_fname)

    dtopo_instant = dtopotools.DTopography(dtopo_instant_fname, dtopo_type=3)


    fig,(ax0,ax1) = plt.subplots(ncols=2,nrows=1,figsize=(12,6))

    X = dtopo_instant.X; Y = dtopo_instant.Y; dZ_at_t = dtopo_instant.dZ_at_t
    dz_max = dtopo_instant.dZ.max()
    tfinal = dtopo_instant.times[-1] + 1  # 1 second after final dZ
    dtopotools.plot_dZ_colors(X,Y,dZ_at_t(tfinal),axes=ax0,
                              cmax_dZ = dz_max,
                              dZ_interval=200, add_colorbar=True);
    ax0.set_title('Seafloor deformation (static seismic)');
    ax0.set_ylim(40,50)
    ax0.set_xlim(-128.5,-122)

    X = dtopo.X; Y = dtopo.Y; dZ_at_t = dtopo.dZ_at_t
    #dz_max = dtopo.dZ.max() # use same color scale as in other figure
    tfinal = dtopo.times[-1] + 1  # 1 second after final dZ
    dtopotools.plot_dZ_colors(X,Y,dZ_at_t(tfinal),axes=ax1,
                              cmax_dZ = dz_max,
                              dZ_interval=200, add_colorbar=True);
    ax1.set_title('Seafloor deformation (static Okada)');
    ax1.set_ylim(40,50)
    ax1.set_xlim(-128.5,-122)

    fname = 'compare_seismic_instant_%s.png' % fault.event
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)



def make_all_cases_okada():
    """
    Output: *caselist*, a list of cases to be run.
    Each case should be dictionary of any parameters needed to set up an
    individual case.  These will be used by run_one_case.

    """

    all_models = \
        ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
         'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']

    #models = all_models
    models = all_models[:3]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

    # Test on one event:
    #events = ['buried-random-NOSUBmur13_deep']
    events = ['buried-random-mur13_deep']

    # Create a list of the cases to be run:
    caselist = []

    for k,event in enumerate(events):
        #event1 = event.replace('buried-','')
        #event_jey = event1.replace('-','_')  # switch to Jey's notation

        fault = set_slip(fault0, event)
        case = {'num':k, 'fault':fault}

        print("+++ in caselist, id(case['fault']) = %i, Mw = %.2f" \
                % (id(case['fault']), case['fault'].Mw()))

        caselist.append(case)

    return caselist


def run_one_case_okada(case):
    """
    Input *case* should be a dictionary with any parameters needed to set up
    and run a specific case.
    """

    num = case['num']
    fault = case['fault']
    event = fault.event
    print('Now running case %i with event %s' % (num,event))
    print('+++ in run_one_case, id(fault) = %i, Mw = %.2f' \
            %  (id(fault),fault.Mw()))
    print('+++ fault.dtopo_times = ', fault.dtopo_times)

    # This part shows how to redirect stdout so output from any
    # print statements go to a unique file...
    import sys
    import datetime

    redirect_output = False

    message = ""
    stdout_fname = 'case%s_out.txt' % case['num']
    try:
        stdout_file = open(stdout_fname, 'w')
        message = message +  "Python output from this case will go to %s\n" \
                            % stdout_fname
    except:
        print(message)
        raise Exception("Cannot open file %s" % stdout_fname)

    print(message) # constructed first to avoid interleaving prints

    if redirect_output:
        sys_stdout = sys.stdout
        sys.stdout = stdout_file
        # send any errors to the same file:
        sys_stderr = sys.stderr
        sys.stderr = stdout_file

    print('\n==============================')
    timenow = datetime.datetime.today().strftime('%Y-%m-%d at %H:%M:%S')
    print('Working on case %s, started at %s' % (case['num'],timenow))

    p = current_process()
    print('Process running this case: ', p)


    fault = case['fault']

    if test_subset:
        # Test on smaller set of subfaults:

        #testfault = dtopotools.Fault(coordinate_specification='triangular')
        testfault = copy.deepcopy(fault)
        testfault.subfaults = []
        for s in fault.subfaults:
            if 45.4<s.latitude<45.6 and -126<s.longitude<-125:
                testfault.subfaults.append(s)
        print('Created testfault with %i subfaults' % len(testfault.subfaults))
        testfault.event = '%s_subset_test' % event
        fault = testfault


    #print('+++ event = ',event)
    #print('+++ event_jey = ',event_jey)
    print('+++ fault.event = ',fault.event)
    print('+++ fault.Mw = %.2f' % fault.Mw())
    print('There are %i subfaults in this model' % len(fault.subfaults))

    slips = array([s.slip for s in fault.subfaults])
    print('+++ max slip = %.2fm' % slips.max())

    if not dry_run:
        dtopo = make_dtopo(fault)  # Apply Okada model
        plot_slip_final_dtopo(fault, dtopo)

        if (not test_subset) and (rupture_type == 'static'):
            try:
                plot_slip_vs_seismic_instant(fault, dtopo)
            except:
                print('plot_slip_vs_seismic_instant failed')

    timenow = datetime.datetime.today().strftime('%Y-%m-%d at %H:%M:%S')
    print('Done with case %s at %s' % (case['num'],timenow))

    if redirect_output:
        # Reset stdout and stdout:
        sys.stdout = sys_stdout
        sys.stderr = sys_stderr

if __name__ == '__main__':

    print('\n--------------------------')
    if dry_run:
        # just print out settings, no runs...
        print('DRY RUN - settings in OkadaStatic.py')

    caselist = make_all_cases_okada()
    nprocs = 1
    run_many_cases_pool(caselist, nprocs, run_one_case_okada)

    if dry_run:
        print('Set dry_run=False and re-execute to make dtopo files')
        print('--------------------------')

    print("Done... See files caseN_out.txt for python output from each case")
