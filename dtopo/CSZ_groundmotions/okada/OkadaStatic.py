#!/usr/bin/env python
# coding: utf-8

# # Adapted from OkadaTest.ipynb 241023
# 
# Read in a fault model from the coarsened triangulation provided by Jey and apply the Okada model from GeoClaw to obtain the static seafloor deformation at the final time.
# 
# (Could be adapted to compute time-dependendent kinetic rupture based on applying Okada at each specified time to the subfaults that have ruptured up to that time, but that would require `rupture time` and `rise_time` for each subfault, which are not included currently in Jey's coarsened model.)
# 


from pylab import *
from clawpack.geoclaw import dtopotools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from copy import copy
import os,sys


# ## Read in the fault geometry:
# 
# And specify a `dtopotools.Fault` object with this geometry.


datadir = '/Users/rjl/D/JeysCode/Audreys_source_models_interpolation_09182024/'
triangles = loadtxt(datadir+'David_cas_fine_mesh.tri')  # list of triangles, 3 vertices
vertices = loadtxt(datadir+'David_cas_fine_mesh.ned')  # list of vertices, lon-lat of each

vertices[:,1] = vertices[:,1] - 360.  # shift to longitude W
vertices[:,3] = -1000*vertices[:,3]   # convert depth to positive depth in meters


# Create Fault with geometry, no slip yet:

fault0 = dtopotools.Fault(coordinate_specification='triangular')
fault0.subfaults = []

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
    import copy
    fault = copy.copy(fault0)
    
    event_dir = datadir + 'output_files_orig'
    
    if 0:
        # don't need slip in dip and strike directions separately, only magnitude read below
        fname = 'dip_slip_resampled_source_saved_%s.out' % event
        dip_slip = loadtxt(os.path.join(event_dir, fname))
        fname = 'strike_slip_resampled_source_saved_%s.out' % event
        strike_slip = loadtxt(os.path.join(event_dir, fname))
        
    fname = 'mag_slip_resampled_source_saved_%s.out' % event
    mag_slip = loadtxt(os.path.join(event_dir, fname))[:,2]
    
    fname = 'rake_resampled_source_saved_%s.out' % event
    rake = loadtxt(os.path.join(event_dir, fname))[:,2]
    
    for j in range(nsubfaults):
        subfault = fault.subfaults[j]
        subfault.rake = rake[j]
        subfault.slip = mag_slip[j]
        
    print('Created fault with Mw = %.2f' % fault.Mw())
    
    # save new name of event in fault object
    fault.event = 'buried-%s_okada_instant' % event.replace('_','-')
    print('New event name: %s' % fault.event)
    return fault


def make_dtopo(fault, times=[0.]):
    if len(times) <= 2:
        fault.rupture_type = 'static'
    else:
        fault.rupture_type = 'kinematic'

    dx = dy = 30/3600.  # spatial resolution for dtopo file
    #x,y = fault.create_dtopo_xy(dx=dx)  # choose automatically
    # use same x,y as in dtopo files made from ground motions:
    x = arange(-128.5, -122.4999, dx)
    y = arange(40,50.0001,dy)
    
    print('Will create dtopo on arrays of shape %i by %i ...' % (len(x),len(y)))
    dtopo = fault.create_dtopography(x,y,times=times,verbose=100);
    
    fname = '%s.dtt3' % fault.event
    dtopo.write(fname, dtopo_type=3)
    print('Created %s with %s displacement at %i times' \
            % (fname, fault.rupture_type, len(times)))
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


all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']
     
models = all_models
#models = all_models[:3]
events = ['%s-deep' % model for model in models] \
       + ['%s-middle' % model for model in models] \
       + ['%s-shallow' % model for model in models]
       
       

# ## Test on one event:
#events = ['locking_mur13_deep']
for event in events:
    event1 = event.replace('buried-','')
    event_jey = event1.replace('-','_')  # switch to Jey's notation
    
    fault = set_slip(fault0, event_jey)
    print('\n==============================')
    print('+++ event = ',event)
    print('+++ event_jey = ',event_jey)
    print('+++ fault.event = ',fault.event)
    print('There are %i subfaults in this model' % len(fault.subfaults))

    if 0:
        dtopo = make_dtopo(fault, times=[0.])
        plot_slip_final_dtopo(fault, dtopo)
        plot_slip_vs_seismic_instant(fault, dtopo)

if 0:
    # Test on smaller set of subfaults:

    testfault = dtopotools.Fault(coordinate_specification='triangular')
    testfault.subfaults = []
    for s in fault.subfaults:
        if 45.4<s.latitude<45.8 and -126<s.longitude<-125:
            testfault.subfaults.append(s)
    print('Created testfault with %i subfaults' % len(testfault.subfaults))
    testfault.event = 'testfault_okada_instant'
    testdtopo = make_dtopo(testfault, times=[0.])
    plot_slip_final_dtopo(testfault, testdtopo)
