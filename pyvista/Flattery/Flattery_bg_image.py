
import os,sys

# import dev_pyvista so that relative imports work:
CLAW = os.environ['CLAW']
VISCLAW = CLAW + '/visclaw'
sys.path.insert(0, VISCLAW)
import dev_pyvista
sys.path.pop(0)

from dev_pyvista.geoclaw.topo_plots.bg_image import fetch_img

extent = [-124.76, -124.55, 48.29, 48.41]
fname = 'Flattery_img.jpg'
fetch_img(extent, fname, source='QuadtreeTiles', scale=16, figsize=None)

if 0:
    fname = 'Flattery_img_map.jpg'
    fetch_img(extent, fname, source='OSM', scale=None, figsize=None)
