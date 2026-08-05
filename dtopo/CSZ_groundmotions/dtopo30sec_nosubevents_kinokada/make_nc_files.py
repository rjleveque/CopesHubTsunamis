
from clawpack.geoclaw import dtopotools
from clawpack.clawutil.util import fullpath_import
import os
from pathlib import Path

CHT = os.environ['CHT']
CHTuser = os.environ['CHTuser']
CHTtools = fullpath_import(f'{CHTuser}/src/CHTuser/CHTtools.py')

dtopo_dir = Path(f'{CHT}/dtopo/CSZ_groundmotions/dtopo30sec_nosubevents_kinokada/dtopofiles')
dtopo_nc_dir = Path(f'{CHT}/dtopo/CSZ_groundmotions/dtopo30sec_nosubevents_kinokada/dtopofiles_nc')

dtopo_nc_dir.mkdir(parents=True, exist_ok=True)

events = CHTtools.all_events()

print(f'Will create {len(events)} nc files in directory:\n   {dtopo_nc_dir}')

for event in events:
    fname = dtopo_dir / f'{event}.dtt3'
    dtopo = dtopotools.DTopography(fname, dtopo_type=3)

    fname_nc_compressed = dtopo_nc_dir / f'{event}.nc'
    dtopo.write(fname_nc_compressed, dtopo_type=4, compression=True)
    print('Created ', f'{event}.nc')


