# *****************************************************************
# *                        PINOCCHIO  V5.1                        *
# *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
# *****************************************************************
# 
# This code was written by
# Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan, 
# Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
# Copyright (C) 2025
# 
# github: https://github.com/pigimonaco/Pinocchio
# web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Directory settings #

RUNDIR = ./

########## System-specific Configuration ############################

SYSTYPE ?= LeonardoBoost

#####################################################################
#####################################################################
# The following variables should be passed from command line to make
# Executable name
EXEC         ?= PinocchioFinal_GPU.x
# Compilation options

# C compiler
COMPILER_CC  ?= mpicc
COMPILER_CPP ?= mpic++
DEBUG        ?= NO
OMP          ?= YES
GPU          ?= NO
FULL_GPU     ?= YES
SCOREP       ?= NO
NSIGHT       ?= NO
ENERGY_CPU   ?= NO
ENERGY_GPU   ?= NO
####################################################################

########## Feature toggles (preprocessor definitions) ##########
#
# How to use:
# - Edit the OPTIONS += -D<FLAG> lines below to enable/disable features.
# - These options enable/disable code paths at compile time. Choose what you need
#   for your run to reduce binary size and keep hot loops lean.
# - The default profile here matches a typical Past Light-Cone (PLC) + Mass Maps
#   run with snapshots and FITS output.
#
# LPT order (displacement field accuracy)
# -DTWO_LPT                    Include second-order LPT term.
# -DTHREE_LPT                  Include third-order LPT term (adds cost).
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# PLC + Mass Maps pipeline
# -DPLC                        Enable past light-cone geometry and basis vectors.
# -DMASS_MAPS                  Build per-segment HEALPix maps counting PLC entries.
# -DHAVE_CFITSIO               Write FITS files (requires CFITSIO installed).
# -DMASS_MAPS_FILTER_UNCOLLAPSED  Keep only entries with Z_PLC > Z_ACC (requires SNAPSHOT).
# -DMASS_MAPS_THREAD_ACCUM     Use per-thread HEALPix accumulators to reduce atomics.
# -DDISP_HISTOGRAM             Optional: write displacement magnitude histograms (debug/QA).
OPTIONS += -DPLC
#OPTIONS += -DMASS_MAPS
#OPTIONS += -DHAVE_CFITSIO
#OPTIONS += -DMASS_MAPS_FILTER_UNCOLLAPSED  # Z_PLC > Z_ACC (requires SNAPSHOT)
#OPTIONS += -DDISP_HISTOGRAM               # Enable only when you need displacement QA
#OPTIONS += -DMASS_MAPS_THREAD_ACCUM

# Mass-maps diagnostics (compile-time)
# -DMASS_MAPS_CAP_DIAG=1        Enable cap prefix diagnostics (off by default)
#   Tracks max accepted/missed angles and prefix misses; useful for aperture QA.
# OPTIONS += -DMASS_MAPS_CAP_DIAG=1

# Triaxial collapse model (choose one)
# -DELL_CLASSIC                Classic ellipsoidal collapse model (default).
# -DELL_SNG                    Spherical-no-gravity approximation.
# -DTABULATED_CT               Use tabulated collapse times from file.
OPTIONS += -DELL_CLASSIC

# Group building and fragmentation
# -DCLASSIC_FRAGMENTATION       Use classic fragmentation instead of default.

# Output and products
# -DSNAPSHOT                   Populate products[] (z_acc, group_ID, displacements, etc.).
#                              Required by MASS_MAPS_FILTER_UNCOLLAPSED.
# -DLIGHT_OUTPUT               Reduce output footprint.
# -DLONGIDS                    Use 64-bit particle IDs.
# -DDOUBLE_PRECISION_PRODUCTS  Store products in double precision.
OPTIONS += -DSNAPSHOT
OPTIONS += -DLONGIDS

# Beyond ΛCDM and solver options
# Neutrinos and scale-dependent growth:
# -DSCALE_DEPENDENT            Enable k-dependent growth (FFT per segment).
# -DREAD_PK_TABLE              Read external P(k) table (implies scale dependence).
#
# Displacement engine:
# -DRECOMPUTE_DISPLACEMENTS    Maintain prev/curr segment fields (needed for PLC crossings).
# -DREAD_HUBBLE_TABLE          Read H(z) table for precise distances.
# -DMASS_MAPS_MATCH_DIAG       Extra cross-matching diagnostics (debug).
# -DUSE_DISTRIBUTE_ALLTOALL    Use all-to-all variant for forward product redistribution.
#OPTIONS += -DMASS_MAPS_MATCH_DIAG
#OPTIONS += -DREAD_PK_TABLE
OPTIONS += -DUSE_DISTRIBUTE_ALLTOALL

# f(R) gravity:
# -DMOD_GRAV_FR                Enable modified gravity model (f(R)).
# -DFR0=1.e-8                  Example FR parameter value.
OPTIONS += -DRECOMPUTE_DISPLACEMENTS
#OPTIONS += -DREAD_HUBBLE_TABLE


# Other options
# -DWHITENOISE                  White-noise ICs.
# -DNORADIATION                 Neglect radiation component.
OPTIONS += -DNORADIATION


# -----------------------------------------------
# - This option impacts CPU code only
# - It is used by default using the GPU with OMP
# OPTIONS += -DCUSTOM_INTERPOLATION

#########################################################################################
#########################################################################################

ifeq ($(SCOREP)$(ENERGY_CPU),YESYES)
$(error "Score-P and PMT are mutually exclusive")
endif

ifeq ($(SCOREP)$(ENERGY_GPU),YESYES)
$(error "Score-P and PMT are mutually exclusive")
endif

################### SCOREP configuration ############################
ifeq ($(SCOREP), YES)
SCOREP_FLAG = scorep --user
OPTIONS += -D_SCOREP
else
SCOREP_FLAG = 
endif
#####################################################################

########## ENERGY OPTIONS ###################################
ENERGY = NO

ifeq ($(ENERGY_CPU), YES)
ENERGY = YES
OPTIONS += -D_ENERGY_RAPL_   # switch-on CPU power measurements through RAPL counters
endif

ifeq ($(ENERGY_GPU), YES)
ENERGY   = YES
OPTIONS += -D_ENERGY_NVIDIA_ # switch-on NVIDIA GPU power measurements
# OPTIONS += -D_ENERGY_AMD_    # switch-on AMD GPU power measurements
endif

ifeq ($(ENERGY), YES)
OPTIONS += -D_ENERGY_PMT_ # use PMT library
COMPILER = $(COMPILER_CPP)
else
COMPILER = $(COMPILER_CC)
endif

#############################################################

######################################################
####################### GPU flag ####################

ifeq ($(GPU), YES)
OPTIONS += -DGPU_OMP
GPU_FLAG = -gpu=ccnative -target=gpu #-gpu=nomanaged
endif

ifeq ($(FULL_GPU), YES)
OPTIONS += -DGPU_OMP_FULL
GPU_FLAG = -gpu=ccnative -target=gpu #-gpu=nomanaged
endif

ifeq ($(NSIGHT), YES)
GPU_FLAG = -gpu=ccnative,lineinfo -target=gpu -gpu=nomanaged
endif

ifeq ($(DEBUG), YES)
# -----------------------------------------------
# - Debug option for multi task MPI run
OPTIONS += -DMPI_ATTACH_DEBUGGER
GPU_FLAG = -gpu=ccnative,debug,lineinfo -target=gpu -gpu=nomanaged 
endif

#######################################################

########## Preprocessor Definitions end ########## 

# OpenMP Configuration #

ifeq ($(OMP), YES)

OMP_FLAG = -fopenmp
OMP_FFTW = -lfftw3_omp -lfftw3_threads

ifeq ($(GPU), YES)
OMP_FLAG = -mp=multicore,gpu
endif

ifeq ($(FULL_GPU), YES)
OMP_FLAG = -mp=multicore,gpu
endif

else
ifeq ($(GPU), YES)
# OMP_FLAG = -mp=gpu
endif

ifeq ($(FULL_GPU), YES)
# OMP_FLAG = -mp=gpu
endif 

OMP_FFTW = -lfftw3_omp
endif # OMP, YES

########################### SYSTYPE CONFIGURATION ###################################

###################### LEONARDO ######################
ifeq ($(SYSTYPE), LeonardoBoost)

CDEBUG      = -g -O0 -Minfo=all -v -Mneginfo $(GPU_FLAG) $(OMP_FLAG)
COPTIMIZED  = -O3 -fast -Minfo=all -v -Mneginfo $(OMP_FLAG)

MPI         = /leonardo/prod/opt/libraries/openmpi/4.1.6/nvhpc--24.3/
MPI_LIBR    = -L$(MPI)/lib -lmpi
MPI_INCL    = -I$(MPI)/include

LIBDIR      = /leonardo_scratch/fast/CNHPC_1498509/lib/nvhpc-23.11

FFTW        = $(LIBDIR)/fftw/fftw-3.3.10
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include

HEFFTE_DIR  = /leonardo/pub/userexternal/dgoz0000/libraries
HEFF_LIBR   = -L$(HEFFTE_DIR)/heffte_nvidia/lib64 -lheffte
HEFF_INCL   = -I$(HEFFTE_DIR)/heffte_nvidia/include

GSL         = $(LIBDIR)/gsl/gsl-2.7.1
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include

PMT         = /leonardo_scratch/fast/CNHPC_1498509/lib/pmt/local
PMT_LIBR    = -L$(PMT)/lib64 -lpmt
PMT_INCL    = -I$(PMT)/include

HEALPIX_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/Healpix_3.83/
CFITSIO_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio


endif

ifeq ($(SYSTYPE), LeonardoDCGP)

CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)

FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW_INC)

GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)

PFFT_ROOT   = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/
PFFT_LIBR   = -L$(PFFT_ROOT)/lib -lpfft
PFFT_INCL   = -I$(PFFT_ROOT)/include

HEALPIX_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/Healpix_3.83/
CFITSIO_ROOT = /leonardo_scratch/large/userexternal/tbatalha/Pinocchio/dep/

HEALPIX_INCL = -I$(HEALPIX_ROOT)/include
HEALPIX_LIBR = -L$(HEALPIX_ROOT)/lib -lchealpix

CFITSIO_INCL = -I$(CFITSIO_ROOT)/include
CFITSIO_LIBR = -L$(CFITSIO_ROOT)/lib -lcfitsio

endif

########## System-specific Configuration end ########## 

# Include paths for libraries #

INC =  $(HEFF_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL) $(PMT_INCL)

# INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)  $(PMT_INCL) # for DCGP

# Library flags #

###############################################################
# HEALPix / CFITSIO configuration (override from environment)
###############################################################

# Append their include paths
INC += $(HEALPIX_INCL) $(CFITSIO_INCL)

# Append their libraries (placed after math deps)
LIB =  -lm $(HEFF_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(CFITSIO_LIBR) $(HEALPIX_LIBR)

# LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(CFITSIO_LIBR) $(HEALPIX_LIBR) # for DCGP

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Source files #

OBJECTS = fmax.o variables.o  cubic_spline_interpolation.o cubic_spline_interpolation_full_GPU.o \
          initialization.o collapse_times.o collapse_times_GPU.o fmax-heffte.o GenIC.o \
	      ReadParamfile.o allocations.o LPT.o distribute.o \
	      fragment.o build_groups.o write_halos.o write_snapshot.o write_mass_maps.o cosmo.o 

ifeq ($(ENERGY), YES)
OBJECTS += energy_pmt_methods.o
endif

# For DCGP
# OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
# 	ReadParamfile.o allocations.o LPT.o distribute.o \
# 	fragment.o build_groups.o write_halos.o write_snapshot.o write_mass_maps.o cosmo.o

# Main targets and rules

pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(SCOREP_FLAG) $(COMPILER) $(COPTS) $(OPTIONS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)
	ldd $(RUNDIR)$(EXEC)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)	

# Generic rule for compiling source files #

energy_pmt_methods.o: energy/energy_pmt_methods.cpp Makefile energy/energy_pmt_methods.h energy/energy_pmt.h
	# g++ -std=c++17 $(OPTIONS) $(PMT_INCL) -c $< -o $@
	$(COMPILER_CPP) -std=c++17 $(OPTIONS) $(PMT_INCL) -c $< -o $@

%.o: %.c Makefile pinocchio.h def_splines.h
	$(SCOREP_FLAG) $(COMPILER_CC) $(INC) $(COPTS) $(OPTIONS)  -c $<

# Clean target #

clean:
	\rm -f *.o *~ $(EXEC) run_planner

###############################################################
# Documentation
###############################################################
# Generates HTML docs from in-source Doxygen comments.
#
# Usage:
#   make docs            # builds HTML docs under ../docs/build/html
#   make docs-clean      # removes generated docs
#
# Requirements:
#   - doxygen installed and available on PATH
#   - config: ../docs/Doxyfile (provided in repo)

.PHONY: docs docs-clean

DOXYGEN ?= doxygen
DOCDIR  ?= ../docs
DOXYCFG ?= $(DOCDIR)/Doxyfile

docs:
	@command -v $(DOXYGEN) >/dev/null 2>&1 || { \
	  echo "Error: doxygen not found. Install doxygen and re-run 'make docs'."; \
	  exit 1; \
	}
	@echo "Generating documentation with $$( $(DOXYGEN) --version )..."
	@cd $(DOCDIR) && $(DOXYGEN) Doxyfile
	@echo "Docs available at $(DOCDIR)/build/html/index.html"

docs-clean:
	@rm -rf $(DOCDIR)/build
