
make .exe
make data  # with no dtopofile
make output

python make_B0_fgmax.py

rsync -avz GH3s_fgmax1_B0.asc $TACC:/home1/04137/rjl/CopesHubTsunamis/geoclaw_runs/onshore_coarse/GH3s_B0/
