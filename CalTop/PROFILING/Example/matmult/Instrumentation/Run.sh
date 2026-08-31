#!/bin/bash

# ==============================================================
# TAU Trace Runner
# ==============================================================

# OpenMP configuration
export OMP_NUM_THREADS=12
export TAU_OMPT_SUPPORT_LEVEL=full

# TAU configuration
export TAU_TRACE=1

# Disable standard profiling output.
# Set to 1 if both profile.* and trace files are desired.
export TAU_PROFILE=0

# Executable
EXECUTABLE="./bin/matrix_multiply.exe"

# Output SLOG2 file
SLOG2_FILE="matrix.slog2"


echo "========================================"
echo "TAU Trace Configuration"
echo "========================================"
echo "Executable       : ${EXECUTABLE}"
echo "OpenMP threads   : ${OMP_NUM_THREADS}"
echo "TAU_TRACE        : ${TAU_TRACE}"
echo "TAU_PROFILE      : ${TAU_PROFILE}"
echo "========================================"


# --------------------------------------------------------------
# Remove old TAU trace files
# --------------------------------------------------------------

echo "Removing old trace files..."

rm -f tautrace.*
rm -f events.*
rm -f tau.trc
rm -f tau.edf
rm -f "${SLOG2_FILE}"


# --------------------------------------------------------------
# Run application
# --------------------------------------------------------------

echo
echo "Running application..."

${EXECUTABLE}

if [ $? -ne 0 ]; then
    echo "ERROR: Application failed."
    exit 1
fi


# --------------------------------------------------------------
# Merge TAU trace files
# --------------------------------------------------------------

echo
echo "Merging TAU traces..."

tau_treemerge.pl

if [ $? -ne 0 ]; then
    echo "ERROR: tau_treemerge.pl failed."
    exit 1
fi


# --------------------------------------------------------------
# Convert TAU trace to SLOG2
# --------------------------------------------------------------

echo
echo "Converting TAU trace to SLOG2..."

tau2slog2 tau.trc tau.edf -o "${SLOG2_FILE}"

if [ $? -ne 0 ]; then
    echo "ERROR: tau2slog2 failed."
    exit 1
fi


# --------------------------------------------------------------
# Open trace with Jumpshot
# --------------------------------------------------------------

echo
echo "Opening trace with Jumpshot..."

jumpshot "${SLOG2_FILE}"