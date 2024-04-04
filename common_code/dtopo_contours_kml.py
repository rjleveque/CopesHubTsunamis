#!/usr/bin/env python
# coding: utf-8

# # make contour plots of dtopo for kml

from pylab import *
from clawpack.geoclaw import topotools, kmltools, dtopotools
from clawpack.visclaw import colormaps, geoplot, gridtools
import os, sys, glob, zipfile

def make_kml(dtopofile, dtopo_type=3, dZ_interval=1, dZmax=20):
    event = os.path.splitext(dtopofile)[0]
    print('Making contours for event = ',event)
    print('  contour interval = %gm' % dZ_interval)
    
    dtopo = dtopotools.DTopography(dtopofile, dtopo_type=dtopo_type)
    
    Zm = ma.masked_where(dtopo.Y < 90, dtopo.Y)  # all masked
    fig,ax,png_extent,kml_dpi = kmltools.pcolorcells_for_kml(dtopo.X, dtopo.Y, Zm,
                                                     png_filename=None, dpc=2)
    clines = arange(dZ_interval,dZmax,dZ_interval)
    lw = 1.
    ax.contour(dtopo.X, dtopo.Y, dtopo.dZ[-1,:,:], clines, colors='r', 
               linestyles='-', linewidths=lw)
    clines = arange(-dZ_interval,dZmax,-dZ_interval)
    ax.contour(dtopo.X, dtopo.Y, dtopo.dZ[-1,:,:], clines, colors='c', 
               linestyles='-', linewidths=lw)
    
    kml_dir = '%s_kml' % event
    os.system('mkdir -p %s' % kml_dir)
    png_filename = kml_dir + '/%s_contours.png' % event
    plt.savefig(png_filename, transparent=True, dpi=kml_dpi)
    
    fname = os.path.join(kml_dir, event+'.kml')
    #fname = name+'.kml'
    png_files = ['%s_contours.png' % event]
    png_names = [event]
    kmltools.png2kml(png_extent, png_files=png_files, png_names=png_names, 
                     name=event, fname=fname)
    
    savedir = os.getcwd()
    os.chdir(kml_dir)
    files = glob.glob('*.kml') + glob.glob('*.png')
    print('kmz file will include:')
    for file in files:
        print('    %s' % os.path.split(file)[-1])
    
    fname_kmz = '%s_contours.kmz' % event
    with zipfile.ZipFile(fname_kmz, 'w') as zip:
        for file in files:
            zip.write(file) 
        print('Created %s' % os.path.abspath(fname_kmz))
    os.chdir(savedir)

if __name__=='__main__':
    dtopofile = sys.argv[1]
    make_kml(dtopofile)
    
