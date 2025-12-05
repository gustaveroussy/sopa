#!/bin/bash

##############################################################
###   This script has to be executed at the root of the    ###
### repository, on a environment with snakemake installed. ###
###             sh tests/test_toy_snakemake.sh             ###
##############################################################

set -e # Exit immediately if any command fails

cd workflow

# Test all toy config files
CONFIG_FILES=(
  config/toy/baysor_overlaps.yaml
  config/toy/baysor_vizgen.yaml
  config/toy/baysor.yaml
  config/toy/cellpose_baysor.yaml
  config/toy/cellpose_comseg.yaml
  config/toy/cellpose_proseg.yaml
  config/toy/cellpose.yaml
  config/toy/comseg.yaml
  config/toy/proseg.yaml
  config/toy/proseg_3d.yaml
  config/toy/backward_comp_check.yaml
  config/toy/cellpose_gpu.yaml
  config/toy/tissue_cellpose.yaml
)


for config in "${CONFIG_FILES[@]}"; do
    echo "🐍 Running snakemake with config: $config"

    rm -rf tuto.zarr tuto.explorer

    snakemake \
        --config sdata_path=tuto.zarr \
        --configfile=$config \
        --workflow-profile profile/local \
        --cores 1

    echo "✅ Pipeline succeeded for $config"
done

echo "🎉 All pipelines ran successfully."
