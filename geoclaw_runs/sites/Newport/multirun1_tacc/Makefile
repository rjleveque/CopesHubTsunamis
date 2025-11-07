# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# on laptop, set SCRATCH to CURDIR so _output is in run directory
#SCRATCH = $(CURDIR)
# executing either one of the lines below will do this ok on laptop
# provided "projects" or "home1" isn't in your path

# on CU cluster this redirects output to scratch disk:
SCRATCH = $(subst /projects/,/rc_scratch/,$(CURDIR))
# on TACC this redirects output to scratch disk:
#SCRATCH = $(subst /home1/,/scratch/,$(CURDIR))

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use

# Environment variable COPEVIS should be set to 
# top of git repository (perhaps in shared locaation on HPC cluster)

# for laptop:
EXE ?= $(CHT)/geoclaw_runs/xgeoclaw_v5-13-1_ifx
#EXE = xgeoclaw


SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = $(SCRATCH)/_output         # Directory for output
SETPLOT_FILE = setplot.py 
PLOTDIR = $(SCRATCH)/_plots         # Directory for plots


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
  $(CHT)/common_code/b4step2.f90 \
  $(CHT)/common_code/getmaxspeed.f90 \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: run all scrdir only_plots check

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
	$(MAKE) data
	$(MAKE) output
	$(MAKE) plots
	#python process_fgout.py

only_plots:
	$(MAKE) plots

check:
	@echo ===================
	@echo running whichclaw shows the following...
	python $(CLAW)/clawutil/src/python/clawutil/whichclaw.py 
	@echo ===================
	@echo CLAW = $(CLAW)
	@echo OMP_NUM_THREADS = $(OMP_NUM_THREADS)
	@echo RUNEXE = $(RUNEXE)
	@echo EXE = $(EXE)
	@echo FC = $(FC)
	@echo FFLAGS = $(FFLAGS)
	@echo LFLAGS = $(LFLAGS)
	@echo OUTDIR = $(OUTDIR)
	@echo PLOTDIR = $(PLOTDIR)
	@echo ===================
