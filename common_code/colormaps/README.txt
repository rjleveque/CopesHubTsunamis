
Colormaps from Yong Wei, 6/24/24

See $CHT/common_code/noaa_colormaps.py
    plot_fgmax_noaa_cmap.py


CHT = os.environ['CHT']   # assuming environment variable set
sys.path.insert(0, os.path.join(CHT,'common_code'))
import noaa_colormaps
cmap_noaa_def, cmap_noaa_max = noaa_colormaps.load()

--------------

from scipy.io import loadmat

mat = loadmat('maxcmap.mat')
maxcmap = mat['maxcmap']   # shape (256,3), columns [R,G,B]

mat = loadmat('defcmap.mat')
defcmap = mat['defcmap']   # shape (256,3), columns [R,G,B]

