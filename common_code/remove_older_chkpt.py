
import pathlib,datetime,glob,contextlib
import numpy as np

dry_run = True

thisdir = pathlib.Path('.')
subdirs = [entry.name for entry in thisdir.iterdir() if entry.is_dir()]
#print('subdirs: ',subdirs)

for subdir in subdirs:

    chkfiles = glob.glob(f'{subdir}/fort.chk*')
    print('in subdir: ',subdir)
    #print('   chkfiles = ',chkfiles)
    if len(chkfiles) > 1:
        fpaths = [pathlib.Path(f) for f in chkfiles]
        timestamps = np.zeros(len(fpaths))
        for i,fpath in enumerate(fpaths):
            timestamp = fpath.stat().st_mtime
            print(f'    {fpath}: {timestamp}')
            timestamps[i] = timestamp
        imax = np.argmax(timestamps)
        for i,fpath in enumerate(fpaths):
            if i != imax:
                print(f'  Will delete {fpath}')
                if not dry_run:
                    fpath.unlink()
                


