
import os,sys

all_models = \
    ['buried-locking-mur13', 'buried-locking-skl16', 'buried-locking-str10',
     'buried-random-mur13',  'buried-random-skl16',  'buried-random-str10']


if 1:

    models = all_models
    events = ['%s-deep' % model for model in models] \
           + ['%s-middle' % model for model in models] \

if 1:
    events = ['buried-locking-str10-middle']

for event in events:
    #srcdir = '_output_%s' % event
    #destdir = 'all_gauges/%s' % event
    srcdir = 'all_gauges/%s' % event
    destdir = 'gnss_data_18events/gnss_data_%s/geoclaw' % event
    os.system('mkdir -p %s' % destdir)
    cmd = 'cp %s/gauge0* %s/' % (srcdir,destdir)
    print(cmd)
    #os.system(cmd)

