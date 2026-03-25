#!/usr/bin/env python
# coding: utf-8

# CopesHubTsunamis/dtopo/CSZ_groundmotions/scenarios_ft_260309/OkadaFromSubsampledSlip.py

"""
Refactored code to handle Frontal Thrust events with 6 frontal faults
labeled A,B,C,D,E,F as well as the megathrust (now called Fault 'M').

To run for a set of events, modify make_cases to specify the set of events
and also specify rupture_type, to make 'static' or 'kinematic'.

"""

from pylab import *
from clawpack.geoclaw import dtopotools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import copy
import os,sys
from multiprocessing import current_process

from clawpack.clawutil.util import fullpath_import
CHTuser = os.environ['CHTuser']
CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set

sys.path.insert(0, os.path.join(root_dir, 'common_code'))
from multip_tools import run_many_cases_pool

dry_run = False  # if True, just print event names
test_subset = False

mu = 30e9  # rigidity = shear modulus (in Pascals)
# Note: the value of mu doesn't matter in Okada solution, only for computing Mw
# Okada depends only on Poisson ratio which is set to 0.25 by default.

# where to find files:
if '/home1/' in os.getcwd():
    # on TACC
    proj_dir = '/corral/projects/NHERI/projects/7f2e74be-d7ca-4e0e-b69a-22c24840b078/CSZ_groundmotions/scenarios_for_tsunami_modeling/'
    geom_dir_M = f'{proj_dir}/scenarios_251229/jey_interp_code_forrandy'
    geom_dir = f'{proj_dir}/source_models/gmsh_faults'
    event_dir = f'{proj_dir}/source_models/jey_interpolated_sources_frontalthrust'
else:
    # RJL laptop
    geom_dir_M = '../scenarios_251229/jey_interp_code_forrandy'
    geom_dir = './gmsh_faults'
    event_dir = './jey_interpolated_sources_frontalthrust'

def load_fault(Nfault, event=None):
    """
    Nfault is a letter, 'M' for megathrust of 'A' thru 'F' for frontal faults

    If event is None, the fault will have subfault geometries set but
    no slip, rake, etc. values.

    This loads one fault and is called repeated from combine_faults defined
    below to create a single set of subfaults for 'M' and the frontal faults.
    """

    fault = dtopotools.Fault(coordinate_specification='triangular')
    fault.subfaults = []

    if Nfault == 'M':
        #geom_dir = '../scenarios_251229/jey_interp_code_forrandy/'
        triangles = loadtxt(f'{geom_dir_M}/JK_cas_fine_mesh_updip_modified.tri')
        vertices = loadtxt(f'{geom_dir_M}/JK_cas_fine_mesh_updip_modified.ned')

        # convert depth to positive depth in meters:
        vertices[:,3] = -1000*vertices[:,3]

    else:
        assert Nfault in 'ABCDEF', f'*** unrecognized Nfault = {Nfault}'
        #geom_dir = './gmsh_faults'
        triangles = loadtxt(f'{geom_dir}/FT{Nfault}.tri')
        vertices = loadtxt(f'{geom_dir}/FT{Nfault}.ned')
        # depth is already positive and in meters

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
        fault.subfaults.append(subfault0)

    print(f'Loaded fault {Nfault} model with {nsubfaults} subfaults')

    if 0:
        # Check that orientation of triangles are all correct:
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


    if event is not None:
        # Also set the slips on the subfaults
        print(f'+++ setting slips for event {event}')

        if Nfault == 'M':
            #event_dir = './jey_interpolated_sources_frontalthrust'

            event_jey = \
                f'{event}_nosubs_scaleMoment_avgabovebelowslab_orig_smoothbelow9900m_correct'
            if 'middle' in event:
                event_jey = event_jey.replace('_correct','_li1cmperyr_correct')

            basename = f'resampled_source_saved_{event_jey}.out'

        else:
            assert Nfault in 'ABCDEF', f'*** unrecognized Nfault = {Nfault}'

            #event_dir = './jey_interpolated_sources_frontalthrust'
            event_dash = event.replace('ft_', 'ft-')
            basename = f'resampled_source_{event_dash}_{Nfault}.out'

        print('+++ in load_fault, loading slips for event = ', event)
        print('+++ event_dir = ',event_dir)
        print('+++ basename = ',basename)

        # don't need slip in dip and strike directions separately, only magnitude

        fname = f'mag_slip_{basename}'
        mag_slip = loadtxt(os.path.join(event_dir, fname))[:,2]

        fname = f'rake_{basename}'
        rake = loadtxt(os.path.join(event_dir, fname))[:,2]

        assert len(rake) == nsubfaults, \
            f'*** For subfault {Nfault}, len(rake) = {len(rake)},' \
            + f'  Expected nsubfaults = {nsubfaults}'

        fname = f'rupt_time_{basename}'
        rupt_time = loadtxt(os.path.join(event_dir, fname))[:,2]
        fname = f'rise_time_{basename}'
        rise_time = loadtxt(os.path.join(event_dir, fname))[:,2]
        rupt_tfinal = (rupt_time + rise_time).max()

        for j in range(nsubfaults):
            subfault = fault.subfaults[j]
            subfault.rake = rake[j]
            subfault.slip = mag_slip[j]
            subfault.mu = mu
            subfault.rupture_time = rupt_time[j]
            subfault.rise_time = rise_time[j]
            subfault.rise_shape = 'quadratic'  # correct?
            subfault.rise_time_starting = None # correct?

        #slips = array([s.slip for s in fault.subfaults])
        print(f'For event {event}, fault {Nfault} has maximum slip ' \
            + f'{abs(mag_slip).max():.2f} and Mw = {fault.Mw():.3f}')

    return fault


def combine_faults(Nfaults='MABCDEF', event=None):
    """
    Make a single fault that includes the subfaults from
        the megathrust (if 'M' in Nfaults)
        plus any frontal faults listed in Nfaults

    If event is None, the fault will have subfault geometries set but
    no slip, rake, etc. values.
    """

    print(f'+++ combining faults {Nfaults}')

    fault = load_fault(Nfaults[0], event)

    for Nfault in Nfaults[1:]:
        next_fault = load_fault(Nfault, event)
        # append all the subfaults of next_fault to those of fault:
        fault.subfaults = fault.subfaults + next_fault.subfaults

    print(f'After combining, fault has {len(fault.subfaults)} subfaults ' \
        + f'with Mw = {fault.Mw():.3f}')

    return fault



def make_dtopo(fault):

    # Grid for dtopo vertical deformations:
    dx = dy = 30/3600.  # spatial resolution for dtopo file
    #x,y = fault0.create_dtopo_xy(dx=dx)  # choose automatically

    # use same x,y as in dtopo files made from ground motions:
    x = arange(-128.5, -122.4999, dx)
    y = arange(40,50.0001,dy)

    print('Will create dtopo on arrays of shape %i by %i ...' % (len(x),len(y)))

    if fault.rupture_type == 'static':
        fault.dtopo_times = [1.]  # only compute static displacement (instant)
    else:
        assert fault.rupture_type in ['kinematic','dynamic'], \
                f'*** Unrecognized rupture_type = {fault.rupture_type}'

        rupt_times = array([s.rupture_time for s in fault.subfaults])
        rise_times = array([s.rise_time for s in fault.subfaults])
        rupt_tfinal = (rupt_times + rise_times).max()
        fault.dtopo_times = arange(0, rupt_tfinal+10, 10)  # every 10 seconds

    # Apply Okada to all subfaults to create deformation dtopo.dZ:
    dtopo = fault.create_dtopography(x,y,times=fault.dtopo_times,verbose=100)

    fname = f'{fault.event}.dtt3'
    dtopo.write(fname, dtopo_type=3)
    print('Created %s with %s displacement at %i times' \
            % (fname, fault.rupture_type, len(dtopo.times)))
    return dtopo


def make_all_cases_okada():
    """
    Output: *caselist*, a list of cases to be run.
    Each case should be dictionary of any parameters needed to set up an
    individual case.  These will be used by run_one_case.

    """

    all_models = \
        ['ft_locking-mur13', 'ft_locking-skl16', 'ft_locking-str10',
         'ft_random-mur13',  'ft_random-skl16',  'ft_random-str10']

    #models = all_models
    models = all_models[:1]
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \
           + ['%s-shallow' % model for model in models]

    # Test on one event:
    #events = ['ft_locking-str10-deep']

    # Create a list of the cases to be run:
    caselist = []

    for k,event in enumerate(events):
        print(f'+++ will combine_faults for event = {event}')
        fault = combine_faults(event=event)
        fault.rupture_type = 'kinematic'

        try:
            fault.event = CHTtools.shortname(event) # e.g. convert to 'FL13D'
        except:
            fault.event = event

        if fault.rupture_type == 'static':
            fault.event = fault.event + '_instant'

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

    # This part shows how to redirect stdout so output from any
    # print statements go to a unique file...
    import sys
    import datetime

    redirect_output = True

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


    print('+++ fault.event = ',fault.event)
    print('+++ fault.Mw = %.3f' % fault.Mw())
    print('There are %i subfaults in this model' % len(fault.subfaults))

    slips = array([s.slip for s in fault.subfaults])
    print('+++ max slip = %.2fm' % slips.max())

    if not dry_run:
        # Apply Okada model:
        dtopo = make_dtopo(fault)

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
        print('DRY RUN - settings in OkadaFromSubsampledSlip.py')

    caselist = make_all_cases_okada()
    print('Will run for events:')
    for case in caselist:
        print(f'    {case['fault'].event}')

    nprocs = 3
    run_many_cases_pool(caselist, nprocs, run_one_case_okada)

    if dry_run:
        print('Set dry_run=False and re-execute to make dtopo files')
        print('--------------------------')

    print("Done... See files caseN_out.txt for python output from each case")
