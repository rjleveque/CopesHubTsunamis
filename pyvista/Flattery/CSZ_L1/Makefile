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

# Environment variable WA_EMD_2020_SHARED should be set to
# top of git repository (perhaps in shared locaation on HPC cluster)

# for laptop:
EXE ?= xgeoclaw

# on CU HPC cluster set EXE as environment variable, so we don't need this:
#EXE = /projects/$(USER)/WA_EMD_Snoho2/cu_cluster/xgeoclaw

SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = $(SCRATCH)/_output         # Directory for output
SETPLOT_FILE = setplot_common.py    # copied from common_python by make all
#SETPLOT_FILE = setplot.py
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
.PHONY: run all scrdir only_plots

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
	cp $(WA_EMD_2020_SHARED)/common_python/setplot_common.py .
	$(MAKE) plots
	python $(WA_EMD_2020_SHARED)/common_python/process_fgmax.py
	python $(WA_EMD_2020_SHARED)/common_python/process_gauges.py

only_plots:
	cp $(WA_EMD_2020_SHARED)/common_python/setplot_common.py .
	$(MAKE) plots
