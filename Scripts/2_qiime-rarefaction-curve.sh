#!/bin/bash

# Activate Qiime2 environment first (line below)
# conda activate qiime2-amplicon-2024.10

# Set up log file
LOGFILE="../Logs/2_qiime_rarefaction-curve.txt"
echo "Starting rarefaction analysis at $(date)" > "$LOGFILE"

# Run qiime diversity alpha-rarefaction (this will take about 1-2 minutes)
{
  echo "Running qiime diversity alpha-rarefaction at $(date)"
  qiime diversity alpha-rarefaction \
    --i-table ../Data/Tables/209766_filtered_feature_table.qza \
    --p-max-depth 3000 \
    --o-visualization ../Data/Qiime2_files/rarefaction_curves.qzv
} >> "$LOGFILE" 2>&1

# Check if the command was successful
if [ $? -eq 0 ]; then
  echo "Rarefaction analysis completed successfully at $(date)" >> "$LOGFILE"
else
  echo "Error during rarefaction analysis at $(date)" >> "$LOGFILE"
fi



### Drag and drop the output qzv file to https://view.qiime2.org/