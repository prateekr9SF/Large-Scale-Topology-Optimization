.PHONY: all clean install uninstall compiler-info

.DEFAULT_GOAL := all

compiler-info:
	@echo "========================================"
	@echo "CalTop compiler configuration"
	@echo "========================================"
	@echo "CC  = $(CC)"
	@$(CC) --version | head -n 1
	@echo
	@echo "CXX = $(CXX)"
	@$(CXX) --version | head -n 1
	@echo
	@echo "FC  = $(FC)"
	@$(FC) --version | head -n 1
	@echo "========================================"

# ============================================================
# Paths for external dependencies
# ============================================================
SPOOLES_PATH = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/Deps/Spooles/src
ARPACK_PATH  = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/Deps/ARPACK
YAML_PREFIX  = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/Deps/yaml-cpp-yaml-cpp-0.6.2/

# Intel MKL
MKL_LIB     = /home/prateek/intel/oneapi/mkl/2025.3/lib
MKL_INCLUDE = /home/prateek/intel/oneapi/mkl/2025.3/include

CMPLR_PATH  = /home/prateek/intel/oneapi/compiler/2025.3/bin


# ============================================================
# TAU Profiling
#
# Normal build:
#   make
#
# TAU profiling build:
#   make PROFILING=1
# ============================================================

PROFILING ?= 0

TAU_ROOT = /home/prateek/Documents/Software/TAU/TAU_install

TAU_MAKEFILE = $(TAU_ROOT)/x86_64/lib/Makefile.tau-caltop-openmp-ompt-pdt-openmp
export TAU_MAKEFILE

TAU_CC  = $(TAU_ROOT)/x86_64/bin/tau_cc.sh
TAU_CXX = $(TAU_ROOT)/x86_64/bin/tau_cxx.sh
TAU_FC  = $(TAU_ROOT)/x86_64/bin/tau_f90.sh

# ============================================================
# Directories
# ============================================================
SRC_DIR     = src
OBJ_DIR     = obj
INCLUDE_DIR = include
BIN_DIR     = bin

# Installation Prefix (Default: Current Directory)
PREFIX ?= $(shell pwd)
BINDIR  = $(PREFIX)/bin

# ============================================================
# Build mode switch
#   MDO=0 : NO preCICE, NO src/Coupling/*, build calTop + calFilt
#   MDO=1 : Enable FSI (preCICE + src/Coupling/*), also build calFSI
# ============================================================
MDO ?= 0

# ============================================================
# Include additional definitions (SCCXF, SCCXC, SCCXCXX, OPT, etc.)
#   Makefile.inc may contain src/Coupling/* sources; we will
#   strip them out below when MDO=0.
# ============================================================
include Makefile.inc

# ============================================================
# Remove all src/Coupling/* sources when MDO=0
# ============================================================
ifeq ($(MDO),0)
  SCCXF   := $(filter-out $(SRC_DIR)/Coupling/%,$(SCCXF))
  SCCXC   := $(filter-out $(SRC_DIR)/Coupling/%,$(SCCXC))
  SCCXCXX := $(filter-out $(SRC_DIR)/Coupling/%,$(SCCXCXX))
endif

# ============================================================
# Compiler and flags
# ============================================================
#CC = cc
#FC = gfortran

CC  = gcc-12
CXX = g++-12
FC  = gfortran-12

LINK_CC = $(CC)
LINK_FC = $(FC)

# ============================================================
# TAU profiling flags
# ============================================================

ifeq ($(PROFILING),1)

  CFLAGS  += -DPROFILING_ON -I$(TAU_ROOT)/include
  #FFLAGS  += -DPROFILING_ON -I$(TAU_ROOT)/include
  
  # PROFILING is 1, use TAU C and FORTRAN wrappers
  LINK_CC = $(TAU_CC)
  LINK_FC = $(TAU_FC)

endif

# Base CFLAGS (no preCICE yet)
CFLAGS = -Wall $(OPT) -fopenmp -pthread \
         -DARCH="Linux" -DARPACK -DPARDISO -DMATRIXSTORAGE -DUSE_MT=1 -DMKL_LP64 \
         -I$(MKL_INCLUDE) -I$(INCLUDE_DIR) 
         #-I$(YAML_PREFIX)/include

# Fortran flags
FFLAGS = -Wall -O3 -fallow-argument-mismatch -std=gnu



# preCICE flags (used only when MDO=1)
PKGCONF_CFLAGS  = $(shell pkg-config --cflags libprecice 2>/dev/null)
PKGCONF_LIBS    = $(shell pkg-config --libs   libprecice 2>/dev/null)

ifeq ($(MDO),1)
  CFLAGS += $(PKGCONF_CFLAGS) -DWITH_PRECICE
else
  CFLAGS += -DNO_PRECICE
endif

# ============================================================
# Source files and objects for main programs
# ============================================================

DRIVER_DIR 	  = $(SRC_DIR)/drivers

SCCXMAIN      = $(DRIVER_DIR)/ccx_2.15.c
SCCXCMAINFILT = $(DRIVER_DIR)/calFilt.c

OCCXMAIN      = $(OBJ_DIR)/drivers/ccx_2.15.o
OCCXMAINFILT  = $(OBJ_DIR)/drivers/calFilt.o

# FSI coupling main (only defined/used when MDO=1)
ifeq ($(MDO),1)
  SCCXMAINFSI = $(SRC_DIR)/Coupling/calFSI.c
  OCCXMAINFSI = $(OBJ_DIR)/Coupling/calFSI.o

  SCCXMAINFSI_ADJ = $(SRC_DIR)/Coupling/calFSI_ADJ.c
  OCCXMAINFSI_ADJ = $(OBJ_DIR)/Coupling/calFSI_ADJ.o
endif


# Extra C sources (you can extend this if needed)
SCCXC_EXTRA = $(SCCXCMAINFILT)
OCCXC_EXTRA = $(SCCXC_EXTRA:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Objects generated from Makefile.inc source lists
OCCXF   = $(SCCXF:$(SRC_DIR)/%.f=$(OBJ_DIR)/%.o)
OCCXC   = $(SCCXC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
OCCXCXX = $(SCCXCXX:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# ============================================================
# Libraries
# ============================================================
LIBS = \
  $(ARPACK_PATH)/libarpack_INTEL.a -Wl,--start-group  \
  $(MKL_LIB)/libmkl_intel_lp64.a \
  $(MKL_LIB)/libmkl_gnu_thread.a \
  $(MKL_LIB)/libmkl_core.a \
  -Wl,--end-group \
  -lpthread -lm -ldl -lstdc++

# YAML-cpp is always linked
LIBS += -L$(YAML_PREFIX)/build -lyaml-cpp

# preCICE libs only when MDO=1
ifeq ($(MDO),1)
  LIBS += $(PKGCONF_LIBS)
endif

# ============================================================
# Default targets (executables)
# ============================================================
ALL_BINS := $(BIN_DIR)/calTop.exe $(BIN_DIR)/calFilt.exe

ifeq ($(MDO),1)
  ALL_BINS += $(BIN_DIR)/calFSI.exe
  ALL_BINS += $(BIN_DIR)/calFSI_ADJ.exe
endif

# Default target
all: $(ALL_BINS)

# Ensure object/bin directories exist
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

# ============================================================
# Object file compilation rules (support subdirectories)
# ============================================================
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) -std=c++11 $(CFLAGS) -c $< -o $@

# ============================================================
# Static library build
#   NOTE: Coupling sources have already been removed from
#   SCCXF/SCCXC/SCCXCXX when MDO=0, so ccx_2.15.a will not
#   contain any src/Coupling/* objects in that case.
# ============================================================
ccx_2.15.a: $(OCCXF) $(OCCXC) $(OCCXCXX) | $(OBJ_DIR)
	ar vr $@ $^

# ============================================================
# Executables
# ============================================================

# calTop
$(BIN_DIR)/calTop.exe: $(OCCXMAIN) ccx_2.15.a | $(BIN_DIR) $(OBJ_DIR)
	./date.pl
	$(CC) $(CFLAGS) -c $(DRIVER_DIR)/ccx_2.15.c -o $(OBJ_DIR)/drivers/ccx_2.15.o
	$(LINK_FC) -Wall -O3 -fopenmp -o $@ $(OCCXMAIN) ccx_2.15.a $(LIBS)
	@echo "Build calTop.exe successfully!"

# calFilt
$(BIN_DIR)/calFilt.exe: $(OCCXMAINFILT) ccx_2.15.a | $(BIN_DIR) $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $(DRIVER_DIR)/calFilt.c -o $(OBJ_DIR)/calFilt.o
	$(LINK_FC) -Wall -O3 -fopenmp -o $@ $(OCCXMAINFILT) ccx_2.15.a $(LIBS)
	@echo "Build calFilt.exe successfully!"

# calFSI (only meaningful when MDO=1)
ifeq ($(MDO),1)
$(BIN_DIR)/calFSI.exe: $(OCCXMAINFSI) ccx_2.15.a | $(BIN_DIR) $(OBJ_DIR)
	$(LINK_FC) -Wall -O3 -fopenmp -o $@ $(OCCXMAINFSI) ccx_2.15.a $(LIBS)
	@echo "Build calFSI.exe successfully!"

# calAdj (only meaningful when MDO=1)
$(BIN_DIR)/calFSI_ADJ.exe: $(OCCXMAINFSI_ADJ) ccx_2.15.a | $(BIN_DIR) $(OBJ_DIR)
	$(LINK_FC) -Wall -O3 -fopenmp -o $@ $(OCCXMAINFSI_ADJ) ccx_2.15.a $(LIBS)
	@echo "Build calFSI_ADJ.exe successfully!"
endif

# ============================================================
# Install / Uninstall
# ============================================================
install: $(ALL_BINS)
	@echo "Installing executables to $(BINDIR)"
	mkdir -p $(BINDIR)
	cp $(BIN_DIR)/calTop.exe  $(BINDIR)/
	cp $(BIN_DIR)/calFilt.exe $(BINDIR)/
ifeq ($(MDO),1)
	cp $(BIN_DIR)/calFSI.exe     $(BINDIR)/
	cp $(BIN_DIR)/calFSI_ADJ.exe $(BINDIR)/
endif
	chmod +x $(BINDIR)/calTop.exe
	chmod +x $(BINDIR)/calFilt.exe
ifeq ($(MDO),1)
	chmod +x $(BINDIR)/calFSI.exe
	chmod +x $(BINDIR)/calFSI_ADJ.exe
endif

uninstall:
	@echo "Removing executables from $(BINDIR)"
	rm -f $(BINDIR)/calTop.exe
	rm -f $(BINDIR)/calFilt.exe
	rm -f $(BINDIR)/calFSI.exe
	rm -f $(BINDIR)/calFSI_ADJ.exe

# ============================================================
# Clean
# ============================================================
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) ccx_2.15.a