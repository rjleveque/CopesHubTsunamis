# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# on UW hyak cluster this redirects output to scratch disk:
SCRATCH = $(subst /mmfs1/home/,/gscratch/tsunami/,$(CURDIR))
RUNDIR = $(CURDIR)

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use

# Environment variable OSVES should be set to 
# top of git repository (perhaps in shared locaation on HPC cluster)

# for laptop:
#EXE ?= $(WP)/xgeoclaw
EXE ?= $(CHT)/Loyce_compliance2028_runs/scenarios_251229/xgeoclaw


SETRUN_FILE = setrun.py
OUTDIR = $(SCRATCH)/_output         # Directory for output
SETPLOT_FILE = setplot.py 
PLOTDIR = $(SCRATCH)/_plots         # Directory for plots


# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
#FFLAGS += -DNETCDF -lnetcdf -I$(NETCDF4_DIR)/include -L$(NETCDF4_DIR)/lib
#LFLAGS += $(FFLAGS) -lnetcdff


# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
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

# Construct the topography data
.PHONY: run all scrdir all_plots

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
	$(MAKE) data > output_makedata.txt
	$(MAKE) output > output.txt

all: 
	$(MAKE) run
	$(MAKE) all_plots

all_plots:
	@echo OUTDIR is $(OUTDIR)
	@echo PLOTDIR is $(PLOTDIR)
	@echo CURDIR is $(CURDIR)
	$(MAKE) plots > output_makeplots.txt
