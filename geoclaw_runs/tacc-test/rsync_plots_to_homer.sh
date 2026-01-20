# Script to rsync the geoclaw_plots directory on the scratch disk
# coupled to this run directory to ptha@homer for viewing on the web.
# Requires ssh keys set up on homer ptha account for user@tacc

# give it an ID that will be used in the directory name on homer:
ID=2026-01-20


#--------------------------------------------
# You shouldn't need to change the rest...
THISDIR=$PWD
RELDIR=${THISDIR//$HOME/}
echo RELDIR = $RELDIR
SCRDIR=${THISDIR//$HOME/$SCRATCH}
echo SCRDIR = $SCRDIR
PLOTDIR=$SCRDIR/geoclaw_plots

REMOTE=public_html${RELDIR}_fromtacc_${ID}
 
echo 
echo -------------------
echo 
echo Will try to copy
echo "  " ${PLOTDIR}
echo to ptha@homer.u.washington.edu:${REMOTE}
echo 
echo This may overwrite any _plots* directories already at:
echo "  "  http://depts.washington.edu/ptha${RELDIR}_fromtacc_${ID}
echo 
echo -n "Is this ok (y/n)? "
read answer
 
if [ "$answer" != "${answer#[Yy]}" ] ;then

    chmod -R og+rX ${PLOTDIR}
    
    rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" \
        ${PLOTDIR} \
        ptha@homer.u.washington.edu:${REMOTE}/
     
    echo 
    echo copied local directory:
    echo "     ${PLOTDIR}"
    echo to homer.  To view plots, open
    echo "     http://depts.washington.edu/ptha${RELDIR}_fromtacc_${ID}/"
     
else
    echo Aborting
fi
 

