# Script to copy a _plots directory to homer for viewing on the web

# fix this path...
RELDIR=CopesHubTsunamis/geoclaw_runs/sites/Newport/multirun2_tacc
PLOTDIR=geoclaw_plots/_plots_BL13D

# give it an ID that will be used in the directory name on homer:
ID=251121a

#--------------------------------------------
# You shouldn't need to change the rest...
 
FULLDIR=/scratch/04137/${USER}/${RELDIR}/${PLOTDIR}
REMOTE=public_html/${RELDIR}/${PLOTDIR}_tacc_${ID}
URL=http://depts.washington.edu/ptha

chmod -R og+rX ${FULLDIR}
 
echo 
echo -------------------
echo    rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" 
echo    ${FULLDIR} 
echo    ptha@homer.u.washington.edu:${REMOTE}
echo -------------------
echo 
echo Will try to copy
echo "  " ${FULLDIR}
echo to ptha@homer.u.washington.edu:${REMOTE}
echo 
echo This will overwrite any _plots already at:
echo "  "  $URL/${RELDIR}/${PLOTDIR}_tacc_${ID}
echo 
echo -n "Is this ok (y/n)? "
read answer
 
if [ "$answer" != "${answer#[Yy]}" ] ;then

    rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" \
        ${FULLDIR} \
        ptha@homer.u.washington.edu:${REMOTE}
     
    echo 
    echo copied local directory:
    echo      ${PLOTDIR}
    echo to homer.  To view plots, open
    echo      ${URL}/${RELDIR}/${PLOTDIR}_tacc_${ID}
     
else
    echo Aborting
fi
 

