#!/bin/bash -e

# Set up an environment to run tests under GitHub Actions
# (see toplevel .github/workflows/build.yml)

IMP_CONDA="imp"

eval "$(conda shell.bash hook)"
conda activate test
conda install --yes -q matplotlib ${IMP_CONDA} pyyaml
