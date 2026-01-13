"""
run in a directory after doing 'make data' so *.data files are
up to date.  Creates a directory B0tempdir to create other files
that are not needed for our purpose, and then it deletes this directory
"""

import os,sys

CHT = os.environ['CHT']

tempdir = 'B0tempdir'
os.system(f'mkdir -p {tempdir}')
os.system(f'cp *.data {tempdir}')

print(f'Running code in temporary directory {tempdir}...')
os.chdir(tempdir)

try:
    os.system(f'{CHT}/common_code/xgeo_computeB')
except:
    print(f'Problem running executable {CHT}/common_code/xgeo_computeB')

os.chdir('..')

print(f'Removing temporary directory {tempdir}')
os.system(f'rm -rf {tempdir}')

