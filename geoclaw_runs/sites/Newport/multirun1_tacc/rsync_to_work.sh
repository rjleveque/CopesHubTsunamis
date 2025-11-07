
RELDIR=CopesHubTsunamis/geoclaw_runs/sites/Newport/multirun1_tacc/geoclaw_plots
PLOTDIR=$SCRATCH/$RELDIR/_plots_buried-locking-mur13-deep
WORKDIR=$WORK/$RELDIR

echo will copy $PLOTDIR
echo      to $WORKDIR

mkdir -p $WORKDIR
rsync -avz $PLOTDIR $WORKDIR/


#rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" \
#        ${PLOTDIR} \
#        ${WORK}/${REMOTE}/
