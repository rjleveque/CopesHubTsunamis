from pylab import *
import os
from scipy.io import loadmat
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from clawpack.visclaw import colormaps,geoplot

# top level directory for this project:
CHT = os.environ['CHT']   # assuming environment variable set
colordir = os.path.join(CHT,'common_code/colormaps')
#sys.path.insert(0, os.path.join(CHT,'common_code'))

def load():
    # load matlab .mat files and convert to RGB arrays:
    mat = loadmat(os.path.join(colordir,'maxcmap.mat'))
    maxcmap = mat['maxcmap']   # shape (256,3)

    mat = loadmat(os.path.join(colordir,'defcmap.mat'))
    defcmap = mat['defcmap']   # shape (256,3)

    #Red = defcmap[:,0]
    #Green = defcmap[:,1]
    #Blue = defcmap[:,2]
    #defcolors = list(zip(Red,Green,Blue))

    #cmap = LinearSegmentedColormap('mycmap',defcolors,N=256)
    #cmap_def = cmap.from_list('mycmap',defcmap)
    #cmap_max = cmap.from_list('mycmap',maxcmap)

    cmap_def = ListedColormap(defcmap)
    cmap_max = ListedColormap(maxcmap)

    return cmap_def, cmap_max

def show(cmap_def, cmap_max):
    figure(301)
    colormaps.showcolors(cmap_def)
    title('(x+y)/2 using cmap_def')
    figure(302)
    colormaps.showcolors(cmap_max)
    title('(x+y)/2 using cmap_max')

