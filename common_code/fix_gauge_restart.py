
import clawpack.pyclaw.gauges as gauges
import copy
import numpy as np
import os, sys, glob


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


def write_gauge(gauge, path=None, format="%+.6e"):
    r"""Write the data from this gauge to a file in `path`

    :Input:
     - *path* (path) Path to write the gauge file to.  Defaults to
       `path = os.getcwd()`.
     - *format* (str) Format string used for the field values.

    :Output:
     None

     Modified from pyclaw.GaugeSolution version to write GeoClaw header, e.g.:
        # gauge_id= 32412 location=( -0.8639200000E+02 -0.1797500000E+02 ) num_var=  4
        # Stationary gauge
        # level, time, q[  1  2  3], eta, aux[]
        # file format ascii, time series follow in this file

    Note: assumes stationary gauge, no aux values written
    Writes time series in .txt file even if original was binary.

    """

    if path is None:
        path = os.getcwd()

    if not gauge.is_valid():
        raise ValueError("Gauge is not initialized properly.")

    gauge_file_name = "gauge%s.txt" % str(gauge.id).zfill(5)
    with open(os.path.join(path, gauge_file_name), "w") as gauge_file:

        gauge_file.write("# gauge_id= %s location=( %s %s ) num_eqn= %s\n" %
             (gauge.id, gauge.location[0], gauge.location[1], gauge.q.shape[0]))
        gauge_file.write("# Stationary gauge\n")
        gauge_file.write("# level, time, q[  1  2  3], eta, aux[]\n")
        gauge_file.write("# file format ascii, time series follow in this file\n")

        # print(gauge.q.shape)
        #import pdb; pdb.set_trace()

        for i in range(gauge.q.shape[1]):
            gauge_file.write("%02i %+.6e " % (gauge.level[i], gauge.t[i]))
            gauge_file.write(" ".join([format % value for value in gauge.q[:, i]]))
            gauge_file.write("\n")

    print(f'Created {path}/{gauge_file_name}')

def fix_all_gauges(outdir='_output', backup=True):

    gauge_files = glob.glob(f'{outdir}/gauge*.txt')
    for gauge_file in gauge_files:
        try:
            gaugeno = int(gauge_file[-9:-4])
        except:
            print(f'+++ ignoring {gauge_file}')
            continue

        orig_gauge_file = gauge_file.replace('.txt','_orig.txt')

        if backup:
            os.system(f'cp {gauge_file} {orig_gauge_file}')
            print(f'Copied {gauge_file} to {orig_gauge_file}')
        gauge = load_gauge_remove_overlap(gaugeno, outdir)
        write_gauge(gauge, path=outdir)

if __name__=='__main__':

    if len(sys.argv) > 1:
        outdir = sys.argv[1]
    else:
        outdir = '_output'

    print('Will overwrite all gauge*.txt files in ', outdir)
    ans = input('OK ? ')
    if ans.strip().lower() == 'y':
        fix_all_gauges(outdir, backup=True)
    else:
        print('Aborting')
