
make .exe
make data  # with no dtopofile
make output

python make_B0_fgmax.py

rsync -avz Newport-Seaside_fgmax1_B0.asc $TACC:/home1/04137/rjl/CopesHubTsunamis/geoclaw_runs/onshore_coarse/Newport-Seaside_B0/
