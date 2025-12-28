from __future__ import print_function

import sys
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')


from pylab import *
import os
from clawpack.geoclaw import topotools, dtopotools

#from inter1d_latitudes_factors import sourcefactors  # use const factor

#buriedpath = 'dtopofiles/buried-locking-str10-deep.dtt3'   #on laptop
#buriedpath = 'dtopofiles/buried-random-str10-deep.dtt3'    #on laptop
WGSpath = 'WGS-DOGAMI-Aexp-B-D-UpS-L-U-2p5.dtt3'            #on laptop

WGS = dtopotools.DTopography()
WGS.read(WGSpath, dtopo_type=3)

WGSnew = dtopotools.DTopography()  # for new dtopo
WGSnew.read(WGSpath, dtopo_type=3)

print('WGS.dZ.shape = ',WGS.dZ.shape)
print('Times in dZ array: ', WGS.times)

assert WGS.dZ.shape[0]==1,      '*** WGS.dZ has different time frames than expected'

dZorig = WGS.dZ
dZorig_lasttime = dZorig[-1,:,:]
print('Original dZ values range from %.2f to %.2f (over all times)' % (dZorig.min(),dZorig.max()))
print('Original dZ values range from %.2f to %.2f (at the last time)'  % (dZorig_lasttime.min(),dZorig_lasttime.max()))

#Plot the original dtopo at time 10s:
figure(figsize=(9,6))
ax1 = subplot(121)
WGS.plot_dZ_colors(t=10, axes=ax1, cmax_dZ=18., dZ_interval=18.)
title('Original WGS \n (At time 10s)')

# increase uplift:
factor_uplift = 1.03
dZnew = where(dZorig>0, dZorig*factor_uplift, dZorig)

# modify subsidence if desired:
factor_subsidence = 1.00
dZnew = where(dZnew<0, dZnew*factor_subsidence, dZnew)

WGSnew.dZ = dZnew

if factor_subsidence==1.:
    fname = 'WGS-DOGAMI-Aexp-B-D-UpS-L-U-2p5_x%.3f.dtt3' % factor_uplift
else:
    fname = 'WGS-DOGAMI-Aexp-B-D-UpS-L-U-2p5_x%.3f_sb_%.3f.dtt3' \
            % (factor_uplift, factor_subsidence)

print('New dZ values range from %.2f to %.2f' % (dZnew.min(),dZnew.max()))
dZnew_lasttime = dZnew[-1,:,:]
print('New dZ values range from %.2f to %.2f (for last time)'  % (dZnew_lasttime.min(),dZnew_lasttime.max()))

WGSnew.write(fname, dtopo_type=3)
print('Created ', fname)

#Plot the new dtopo at time 10s:
fnprint = fname + '\n (at time 10s)'
ax2 = subplot(122)
WGSnew.plot_dZ_colors(t=10, axes=ax2, cmax_dZ=18., dZ_interval = 18.)
title(fnprint)
fname_png = os.path.splitext(fname)[0] + '.png'
savefig(fname_png)
print('Created ', fname_png)

#Plot differences in dtopo at time t=10s:
figure(figsize=(9,6))
dWGS = dtopotools.DTopography() # just for plotting diffs
dWGS.read(WGSpath, dtopo_type=3)
dWGS.dZ = dZnew - dZorig
dWGS.plot_dZ_colors(t=10, cmax_dZ=2., dZ_interval=1.)
title('Change in dZ %s' % fnprint)
fname_png = 'difference_to_WGS-DOGAMI-Aexp-B-D-UpS-L-U-2p5.png'
savefig(fname_png)
print('Created ', fname_png)
#########

