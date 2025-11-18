# Script to copy file or directory to the MyData/stampede3_results directory
# on DesignSafe for processing on the Jupyter Hub.
# Does not copy fort.* files (frame output or checkpoint files)

# fix this path...
RELDIR=CopesHubTsunamis/geoclaw_runs/sites/GraysHarborBridges/WesportBridges_multirun_tacc/geoclaw_outputs
TOCOPY=_output_BL10D

# give it an ID that will be used in the directory name on homer:
#ID=251117a

#--------------------------------------------
# You shouldn't need to change the rest...
 
FULLDIR=/scratch/04137/${USER}/${RELDIR}
XFERDIR=stampede3_results
REMOTE=${USER}@designsafe.data.tacc.utexas.edu:/data/designsafe/mydata/${USER}/${XFERDIR}/${RELDIR}
 
echo 
echo -------------------
echo 
echo Will try to copy
echo "  " ${FULLDIR}/${TOCOPY}
echo to MyData/${XFERDIR}/
echo 
echo This will overwrite any existing data
echo 
echo -n "Is this ok (y/n)? "
read answer
 
if [ "$answer" != "${answer#[Yy]}" ] ;then

    rsync -avz --rsync-path="mkdir -p ${REMOTE} && rsync" \
	--exclude 'fort.*' \
        ${FULLDIR}/${TOCOPY} ${REMOTE}/
     
    echo 
    echo copied local directory:
     
else
    echo Aborting
fi
 

