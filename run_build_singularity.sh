#!/bin/bash

# Export temporary build paths to local disk
export SINGULARITY_TMPDIR=/tmp
export SINGULARITY_CACHEDIR=/tmp

TARGET_DIR="/archive/InternalMedicine/Chung_lab/shared/sc/sig_image"
mkdir -p "${TARGET_DIR}"

# Build image
singularity build --fakeroot "${TARGET_DIR}/pipeline_tools.sif" singularity.def

#TARGET_DIR="/archive/InternalMedicine/Chung_lab/shared/sc/sig_image"
#TARGET_DIR="./"
#mkdir -p "${TARGET_DIR}"
#singularity build --fakeroot "${TARGET_DIR}/pipeline_tools.sif" singularity.def
