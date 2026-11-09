
Where to find dtopo files:

On laptop:
    Subdirectories of $CHT/dtopo/CSZ_groundmotions

On stampede3: Subdirectories of
    /work2/04137/rjl/CHTshare/CopesHubTsunamis/dtopo/CSZ_groundmotions

On DesignSafe:
    https://www.designsafe-ci.org/data/browser/projects/PRJ-6005/
        Working Directory --> CSZ_groundmotions
    or accessible from login node on stampede3 as:
    /corral/projects/NHERI/projects/7f2e74be-d7ca-4e0e-b69a-22c24840b078/CSZ_groundmotions/

# Original seismic sources (F events updated 11/6/26):
dtopo30sec_seismic/
    dtopofiles/
        BL10D.dtt3 # etc

# From KinOkada on modified slips with no subevents:
dtopo30sec_nosubevents_kinokada/
    dtopofiles/
        BL10D.dtt3 # etc 
    dtopofiles_instant/
        BL10D_instant.dtt3 # etc 
    dtopofiles_nc/
        BL10D.nc # etc 
