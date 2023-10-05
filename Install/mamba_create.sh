#!/bin/bash
#
# ---- user parameters ----
# set path for mamba binary and Pyticles environment path
mamba_exe="/home2/datahome/jcollin/.local/bin/micromamba"
mamba_root_prefix="/home2/datawork/jcollin/conda"
pyticles_env="pyticles-3.11.yml"

# ---- end of user parameters -----
export MAMBA_EXE=$mamba_exe
export MAMBA_ROOT_PREFIX=$mamba_root_prefix

$mamba_exe env create -f $pyticles_env 

