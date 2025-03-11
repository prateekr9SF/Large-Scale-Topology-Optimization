.PHONY: all clean install uninstall

# Define paths for SPOOLES and ARPACK
SPOOLES_PATH = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/Deps/Spooles/src
ARPACK_PATH = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/Deps/ARPACK

#Directories
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include
BIN_DIR = bin

# Installation Prefix (Default: Current Directory)
PREFIX ?= $(shell pwd)
BINDIR = $(PREFIX)/bin



# Compiler and flags
CFLAGS = -Wall -O3 -I $(SPOOLES_PATH) -I $(INCLUDE_DIR) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O3 -fallow-argument-mismatch

CC = cc
FC = gfortran



# Defauilt target
all: $(BIN_DIR)/calTop.exe

# Ensure object directory exists
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@


# Include additional makefile if necessary
include Makefile.inc

# Source files and objects
SCCXMAIN = $(SRC_DIR)/ccx_2.15.c

OCCXF = $(SCCXF:$(SRC_DIR)/%.f=$(OBJ_DIR)/%.o)
OCCXC = $(SCCXC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
OCCXMAIN = $(OBJ_DIR)/ccx_2.15.o

# Libraries
LIBS = \
	$(SPOOLES_PATH)/spooles.a \
	$(ARPACK_PATH)/libarpack_INTEL.a \
	-lpthread -lm -lc


# Object file compilation rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Build executable in bin/
$(BIN_DIR)/calTop.exe: $(OCCXMAIN) ccx_2.15.a $(LIBS) | $(BIN_DIR) $(OBJ_DIR)
	./date.pl
	$(CC) $(CFLAGS) -c $(SRC_DIR)/ccx_2.15.c -o $(OBJ_DIR)/ccx_2.15.o
	$(FC) -Wall -O3 -o $@ $(OCCXMAIN) ccx_2.15.a $(LIBS)
	@echo "Build successful!"

# Build static library
ccx_2.15.a: $(OCCXF) $(OCCXC) | $(OBJ_DIR)
	ar vr $@ $^


# Install target (default: install in current directory)
install: $(BIN_DIR)/calTop.exe
	@echo "Installing calTop.exe to $(BINDIR)"
	mkdir -p $(BINDIR)
	cp $(BIN_DIR)/calTop.exe $(BINDIR)/
	chmod +x $(BINDIR)/calTop.exe


# Uninstall target
uninstall:
	@echo "Removing calTop.exe from $(BINDIR)"
	rm -f $(BINDIR)/calTop.exe

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/*.* ccx_2.15.a
