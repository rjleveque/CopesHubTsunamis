# Script to copy a _plots directory to homer for viewing on the web

# fix this path...
RELDIR=CopesHubTsunamis/geoclaw_runs/CSZ_groundmotions/X

# give it an ID that will be used in the directory name on homer:
ID=2024-04-13

#--------------------------------------------
# You shouldn't need to change the rest...
 
FULLDIR=/scratch/04137/${USER}/${RELDIR}
PLOTDIR=${FULLDIR}/_plots
REMOTE=public_html/${RELDIR}_tacc_${ID}
 
echo 
echo -------------------
echo 
echo Will try to copy
echo "  " ${PLOTDIR}
echo to ptha@homer.u.washington.edu:${REMOTE}
echo 
echo This will overwrite any _plots already at:
echo "  "  http://depts.washington.edu/ptha/${RELDIR}_hyak_${ID}
echo 
echo -n "Is this ok (y/n)?"
read answer
 
if [ "$answer" != "${answer#[Yy]}" ] ;then

    rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" \
        ${PLOTDIR} \
        ptha@homer.u.washington.edu:${REMOTE}/
     
    echo 
    echo copied local directory:
    echo      ${PLOTDIR}
    echo to homer.  To view plots, open
    echo      http://depts.washington.edu/ptha/${RELDIR}_hyak_${ID}/_plots/_PlotIndex.html
     
else
    echo Aborting
fi
 

