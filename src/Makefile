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
EXEC         ?= hinocchio.x
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

########## Preprocessor Definitions ################################

# Displacement LPT order # 

OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# PLC reconstruction

# OPTIONS += -DPLC

# Dynamics of triaxial collapse #

OPTIONS += -DELL_CLASSIC
# OPTIONS += -DELL_SNG
# OPTIONS += -DTABULATED_CT

# Building groups and fragmentation #

# OPTIONS += -DCLASSIC_FRAGMENTATION

# Output #

# OPTIONS += -DSNAPSHOT
# OPTIONS += -DLIGHT_OUTPUT
# OPTIONS += -DLONGIDS
# OPTIONS += -DDOUBLE_PRECISION_PRODUCTS

# Beyond LambdaCDM models

# these are for neutrino cosmology
# OPTIONS += -DSCALE_DEPENDENT
# OPTIONS += -DREAD_PK_TABLE
# OPTIONS += -DONLY_MATTER_POWER
# OPTIONS += -DRECOMPUTE_DISPLACEMENTS

# add also these for f(R) gravity
# OPTIONS += -DMOD_GRAV_FR
# OPTIONS += -DFR0=1.e-8

# Other options # 

# OPTIONS += -DWHITENOISE
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
GPU_FLAG = 

ifeq ($(GPU), YES)
OPTIONS += -DGPU_OMP
GPU_FLAG = -gpu=ccnative -target=gpu #-gpu=nomanaged
endif

ifeq ($(FULL_GPU), YES)
OPTIONS += -DGPU_OMP_FULL
GPU_FLAG = -gpu=ccnative -target=gpu #-gpu=nomanaged
GPU_FLAG =
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
OMP_FLAG = -mp=multicore # -fopenmp for gcc/clang!
OMP_FFTW = -lfftw3_omp -lfftw3_threads

ifeq ($(GPU), YES)
# OMP_FLAG = -mp=multicore,gpu
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

OMP_FLAG =
OMP_FFTW = -lfftw3_omp
endif # OMP, YES

########################### SYSTYPE CONFIGURATION ###################################

###################### LeonardoBoost ######################

ifeq ($(SYSTYPE),LeonardoBoost)
CDEBUG      = -g -O0 -Minfo=all -v -Mneginfo $(GPU_FLAG) $(OMP_FLAG)
COPTIMIZED  = -O3 -fast -march=native $(GPU_FLAG) $(OMP_FLAG)
LIBDIR      = /leonardo/pub/userexternal/dgoz0000/libraries
FFTW        = $(LIBDIR)/my_fftw_nvidia
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include
GSL         = /leonardo_scratch/fast/CNHPC_1498509/lib/nvhpc-23.11/gsl/gsl-2.7.1
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include
HEFF_LIBR   = -L$(LIBDIR)/heffte_nvidia/lib64 -lheffte
HEFF_INCL   = -I$(LIBDIR)/heffte_nvidia/include
endif

##########################################################################


# ###################### PLEIADI ######################

ifeq ($(SYSTYPE), pleiadi)

COPTIMIZED  = -O3 -g -fopenmp -march=native -mtune=native -ftree-vectorize
FFTW        = /u/dgoz/lib/fftw-3.3.10_pleiadi
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include
PFFT        = /u/dgoz/lib/pfft-1.0.8-alpha_pleiadi
PFFT_LIBR   = -L$(PFFT)/lib -lpfft
PFFT_INCL   = -I$(PFFT)/include
MPI_LIBR    =
MPI_INCL    =
GSL         = /u/dgoz/lib/gsl-2.7.1_pleiadi
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include

endif
######################################################

#MAGELLANUS

######################################################

ifeq ($(SYSTYPE),"magellanus")
CC          =  mpicc 
CDEBUG      =  -Wall -g #$(OMP_FLAG)
#COPTIMIZED  =  -g -O3 -march=native -mp=multicore,gpu -gpu=ccnative -target=gpu #$(OMP_FLAG)
COPTIMIZED  =  -g -O3 -march=native $(OMP_FLAG) $(GPU_FLAG)
#COPTIMIZED  =  -g -O3 -march=native -mp=multicore $(OMP_FLAG)
#COPTIMIZED  =  -g -O0 -march=native
FFTW_LIBR   =  -L/data/NVC_LIB/FFT_lib/fftw-3.3.10/local/lib $(OMP_FFTW)
FFTW_INCL   = -I/data/NVC_LIB/FFT_lib/fftw-3.3.10/local/include
HEFF_LIBR   = -L/data/NVC_LIB/FFT_lib/heffte/local/lib -lheffte
HEFF_INCL   = -I/data/NVC_LIB/FFT_lib/heffte/local/include
PFFT_LIBR   = -L/data/NVC_LIB/FFT_lib/pfft/pfft-1.0.8-alpha/local/lib -lpfft
PFFT_INCL   = -I/data/NVC_LIB/FFT_lib/pfft/pfft-1.0.8-alpha/local/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L/data/NVC_LIB/gsl-2.8/local/lib -lgsl -lgslcblas
GSL_INCL    = -I/data/NVC_LIB/gsl-2.8/local/include
endif

########## System-specific Configuration end ########## 

# Include paths for libraries #

INC =  $(HEFF_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL) $(PMT_INCL)

# Library flags #

LIB = -lm $(HEFF_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(PMT_LIBR)

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Print GPU_FLAG and COPTS for debugging
$(info GPU_FLAG = $(GPU_FLAG))
$(info COPTS    = $(COPTS))


# Source files #

OBJECTS = fmax.o  variables.o  cubic_spline_interpolation.o cubic_spline_interpolation_full_GPU.o cosmo.o initialization.o fmax-heffte.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o collapse_times.o collapse_times_GPU.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o pinocchio.o mpi_debugger.o

ifeq ($(ENERGY), YES)
OBJECTS += energy_pmt_methods.o
endif

# Main targets and rules

$(EXEC): $(OBJECTS) Makefile
	$(SCOREP_FLAG) $(COMPILER) $(COPTS) $(OPTIONS) -o $(RUNDIR)$(EXEC) $(OBJECTS) $(LIB)
	ldd $(RUNDIR)$(EXEC)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)	

# Generic rule for compiling source files #

energy_pmt_methods.o: energy/energy_pmt_methods.cpp Makefile energy/energy_pmt_methods.h energy/energy_pmt.h
	# g++ -std=c++17 $(OPTIONS) $(PMT_INCL) -c $< -o $@
	$(COMPILER_CPP) -std=c++17 $(OPTIONS) $(PMT_INCL) -c $< -o $@

%.o: %.c Makefile pinocchio.h def_splines.h
	$(SCOREP_FLAG) $(COMPILER_CC) $(INC) $(COPTS) $(OPTIONS) -c $<


# Clean target #

clean:
	rm -f *.o *~ $(EXEC) run_planner
