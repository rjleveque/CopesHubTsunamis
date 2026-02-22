# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Environment variable CHT should be set to point to the
# top of git repository, in order to use the same executable in all directories.
CHT ?= $(CURDIR)
#EXE ?= $(CHT)/geoclaw_runs/xgeoclaw-share
EXE = xgeo
# Note ?= since on HPC machine EXE might be set in slurm script instead


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
SETRUN_FILE = setrun_case.py        # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot_case.py      # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------


MODULES = \

SOURCES = \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

.PHONY: all scrdir
	
scrdir: 
	mkdir -p $(SCRATCH)
	@echo Using data in RUNDIR = $(CURDIR)
	@echo Sending output to OUTDIR = $(OUTDIR)
	@echo "OUTDIR on $$(date):" >> OUTDIR.txt
	@echo "$(OUTDIR)" >> OUTDIR.txt
	@echo "RUNDIR on $$(date):" >> $(SCRATCH)/RUNDIR.txt
	@echo "$(CURDIR)" >> $(SCRATCH)/RUNDIR.txt

run: 
	$(MAKE) scrdir
	$(MAKE) data
	$(MAKE) output
	
all:
	$(MAKE) scrdir 
	$(MAKE) .plots
	$(MAKE) .htmls

