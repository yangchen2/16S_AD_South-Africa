# Step 1: Convert HDF5 → JSON BIOM
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all_json.biom \
  --to-json

# Step 2: Convert JSON BIOM → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all.tsv \
  --to-tsv

# You will then need to get rid of the header row and the index name (start with #) before you run the plotting notebooks