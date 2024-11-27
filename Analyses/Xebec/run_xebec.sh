#!/bin/bash -l
# Bash script to run Xebec workflow
# Author: Yang Chen
# Date: 11-27-2024
# Description: Executes Xebec with specified input files and outputs.

# Activate the conda environment
echo "Activating xebec_python39 environment..."
source activate xebec_python39

# Run Xebec command
echo "Running Xebec..."
xebec \
  -ft ../../Data/Tables/209723_reference-hit_nwk-matched.biom \
  -m ../../Data/Metadata/updated_clean_ant_skin_metadata.tab \
  -t ../../Data/Trees/209723_insertion_tree.relabelled.nwk \
  -o Output/

echo "Xebec execution completed."

