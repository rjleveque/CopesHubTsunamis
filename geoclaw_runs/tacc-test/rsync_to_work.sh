#!/bin/bash

# Copy plots from $SCRATCH to $WORK so they can be viewed on 
# the DesignSafe JupyterHub
#
# This assumes job run from this directory on $HOME, e.g.
#   /home1/04137/rjl/CopesHubTsunamis/geoclaw_runs/tacc-test
# created a scratch directory for output and plots with the
# same directory structure on $SCRATCH, e.g.
#   /scratch/04137/rjl/CopesHubTsunamis/geoclaw_runs/tacc-test/geoclaw_plots
# and we want to copy to the same directory structure on $WORK, e.g.
#   /work2/04137/rjl/stampede3/CopesHubTsunamis/geoclaw_runs/tacc-test/geoclaw_plots
#
# Then on DesignSafe, open the plots in 
#  Work/CopesHubTsunamis/geoclaw_runs/tacc-test/geoclaw_plots

THISDIR=$PWD
RELDIR=${THISDIR//$HOME/}
#echo RELDIR = $RELDIR
SCRDIR=${THISDIR//$HOME/$SCRATCH}
#echo SCRDIR = $SCRDIR
PLOTDIR=$SCRDIR/geoclaw_plots
WORKDIR=$WORK/$RELDIR/geoclaw_plots

echo will copy -r $PLOTDIR
echo      to      $WORKDIR

mkdir -p $WORKDIR
rsync -avz $PLOTDIR/ $WORKDIR/

