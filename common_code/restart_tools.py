

def find_last_checkpt(outdir, verbose=True):
    import glob
    import numpy as np
    forttck_files = glob.glob(f'{outdir}/fort.tck*')
    if len(forttck_files) == 0:
        if verbose:
            print('Found no checkpoint file')
        restart_file = None
        return restart_file
        
    tchks = [float(open(f,'r').readline().split()[-1]) for f in forttck_files]
    if verbose:
        print('Found checkpoints at times: ', tchks)

    k = np.argmax(tchks)
    forttck = forttck_files[k]
    restart_file = forttck.replace('tck','chk')
    return restart_file
        

def time(restart_file):
    forttck = restart_file.replace('chk','tck')
    t = float(open(forttck,'r').readline().split()[-1])
    return t

