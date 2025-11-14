############################################################
# BIOM → JSON → TSV CONVERSION PIPELINE
# Converts multiple BIOM tables (skin/nasal, 1%/10% prev-filt)
# into TSV format for downstream plotting.
#
# NOTE:
# After conversion, open the TSV files and remove:
#   - the first header row starting with "#"
#   - the index name row
############################################################


############################################################
# 10% PREVALENCE FILTER – ALL SAMPLES
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all_json.biom \
  -o ../Data/Tables/Count_Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all.tsv \
  --to-tsv


############################################################
# 10% PREVALENCE FILTER – SKIN ONLY
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_skin.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_skin_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_skin_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_skin.tsv \
  --to-tsv


############################################################
# 10% PREVALENCE FILTER – NASAL ONLY
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_nasal.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_nasal_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_nasal_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_nasal.tsv \
  --to-tsv


############################################################
# 1% PREVALENCE FILTER – ALL SAMPLES
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_all.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_all_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_all_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_all.tsv \
  --to-tsv


############################################################
# 1% PREVALENCE FILTER – NASAL ONLY
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_nasal.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_nasal_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_nasal_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_nasal.tsv \
  --to-tsv


############################################################
# 1% PREVALENCE FILTER – SKIN ONLY
############################################################

# HDF5 → JSON
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_skin.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_skin_json.biom \
  --to-json

# JSON → TSV
biom convert \
  -i ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_skin_json.biom \
  -o ../Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-1pct_Genus-ASV_skin.tsv \
  --to-tsv


