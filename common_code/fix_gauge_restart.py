
import clawpack.pyclaw.gauges as gauges
import copy
import numpy as np


def load_gauge_remove_overlap(gaugeno, outdir):
    """
    read in gauge data for one gaugeno and filter it to remove any
    overlapping times that resulted from a restart, keeping the 
    data newer data.

    Assumes there is only one overlapping section to remove.

    returns the filtered data as a gauges.GaugeSolution object
    """

    gauge = gauges.GaugeSolution(gaugeno, outdir)
    gauge_new = copy.copy(gauge)

    t_orig = gauge.t
    tdiff = np.diff(t_orig)
    if tdiff.min() >= 0:
        print(f'gauge {gaugeno} t values are non-decreasing, not changing')
    else:
        inew1 = np.where(tdiff < 0)[0][0]+1  # index where first overlap starts
        tnew1 = t_orig[inew1]   # time new data in overlap starts
        inew2 = len(t_orig) - 1 # final index, assuming one overlap
        tnew2 = t_orig[inew2]   # time new data in overlap ends
        iold2 = inew1           # last index of section to remove
        told2 = t_orig[iold2]
        iold1 = np.where(t_orig >= tnew1)[0][0]  # first index to remove
        told1 = t_orig[iold1]
        print(f'Will remove indices {iold1} to {iold2} from time {told1:.3f} to' \
              f' {told2:.3f}')
        print(f'Remaining times will go from {tnew1} to {tnew2}')

        gauge_new.t = np.hstack((t_orig[:iold1], t_orig[iold2:]))
        gauge_new.q = np.hstack((gauge.q[:,:iold1], gauge.q[:,iold2:]))
        
    #return gauge, gauge_new  # for debugging

    return gauge_new

