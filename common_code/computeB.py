"""
run in a directory after doing 'make data' so *.data files are
up to date.  Creates a directory B0junkdir to create other files
that are not needed for our purpose, and then it deletes this directory
"""

import os,sys

CHT = os.environ['CHT']

junkdir = 'B0junkdir'
os.system(f'mkdir -p {junkdir}')
os.system(f'cp *.data {junkdir}')
os.chdir(junkdir)

os.system(f'{CHT}/common_code/xgeo_computeB')

os.chdir('..')
os.system(f'rm -rf {junkdir}')

