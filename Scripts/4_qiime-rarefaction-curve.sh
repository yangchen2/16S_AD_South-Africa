#!/bin/bash

# Activate Qiime2 environment first (line below)
# conda activate qiime2-amplicon-2024.10

##########################################################################################
# SCRIPT 4: QIIME2 FEATURE TABLE SUMMARY AND RAREFACTION CURVE
##########################################################################################

# Set up log file
LOGFILE="../Logs/4_qiime_rarefaction-curve.log"
echo "Starting rarefaction analysis at $(date)" > "$LOGFILE"

# Summarize the feature table to inspect sample read depths
{
  echo "Running qiime feature-table summarize at $(date)"
  qiime feature-table summarize \
    --i-table ../Data/Tables/Count_Tables/3_209766_feature_table_dedup_prev-filt-1pct.qza \
    --o-visualization ../Data/QC/feature_table_summary.qzv
    
  echo "Feature table summary created: ../Data/QC/feature_table_summary.qzv"
  echo "Use https://view.qiime2.org/ to view per-sample sequencing depths before choosing rarefaction cutoff."
} >> "$LOGFILE" 2>&1


# Run qiime diversity alpha-rarefaction (this will take about 1-2 minutes)
{
  echo "Running qiime diversity alpha-rarefaction at $(date)"
  qiime diversity alpha-rarefaction \
    --i-table ../Data/Tables/Count_Tables/3_209766_feature_table_dedup_prev-filt-1pct.qza \
    --p-max-depth 3000 \
    --o-visualization ../Data/QC/rarefaction_curves.qzv
} >> "$LOGFILE" 2>&1

# Check if the command was successful
if [ $? -eq 0 ]; then
  echo "Rarefaction analysis completed successfully at $(date)" >> "$LOGFILE"
else
  echo "Error during rarefaction analysis at $(date)" >> "$LOGFILE"
fi

### Drag and drop the output qzv file to https://view.qiime2.org/