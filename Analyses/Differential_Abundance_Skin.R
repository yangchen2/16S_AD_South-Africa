###########################################################################
# Set working directory
setwd("/Users/yangchen/PhD/Collaborations/Dube_lab/16S_AD_South-Africa/Analyses")
# Load packages
library(ANCOMBC)
library(tidyverse)
library(biomformat)
library(phyloseq)
cat("ANCOM-BC2 version:", as.character(packageVersion("ANCOMBC")), "\n")
###########################################################################
# Define paths
metadata_path <- '../Metadata/16S_AD_South-Africa_metadata_subset.tsv'
biom_path <- '/Users/yangchen/PhD/Collaborations/Dube_lab/16S_AD_South-Africa/Data/Tables/Count_Tables/6_209766_feature_table_dedup_prev-filt-10pct_Genus-ASV_all.tsv'


output_dir <- "../Data/Differential_Abundance"
dir.create(output_dir, showWarnings = FALSE)
###########################################################################
# Load metadata
metadata <- read_tsv(metadata_path)
metadata <- metadata %>%
  rename(sample_id = `#sample-id`) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", "")) %>%
  column_to_rownames("sample_id")

cat("\nMetadata sample IDs (first 20):\n")
print(head(rownames(metadata), 20))
###########################################################################
# Load OTU table
otu_data <- read_tsv(biom_path, comment = '#') %>% column_to_rownames(colnames(.)[1])

cat("\nOTU table row names (first 20):\n")
print(head(colnames(otu_data), 20))

cat("OTU table loaded:", ncol(otu_data), "samples,", nrow(otu_data), "features\n")
###########################################################################
# Define reusable function
run_ancombc2_by_area <- function(area_name, metadata, otu_data, output_dir) {
  cat("\n==========================\nRunning for:", area_name, "\n==========================\n")
  
  # Filter metadata
  meta_sub <- metadata %>% filter(area == area_name)
  if (nrow(meta_sub) < 10) {
    message("Skipping ", area_name, " — fewer than 10 samples.")
    return(NULL)
  }
  
  cat("Samples in", area_name, ":", nrow(meta_sub), "\n")

  # Align sample IDs
  common_samples <- intersect(rownames(meta_sub), colnames(otu_data))
  meta_sub <- meta_sub[common_samples, ]
  otu_sub <- otu_data[, common_samples]
  cat("After alignment:", nrow(meta_sub), "samples\n")
  
  if (nrow(meta_sub) < 10) {
    message("Skipping ", area_name, " — fewer than 10 aligned samples.")
    return(NULL)
  }
  
  # Check for missing values in covariates
  covariate_cols <- c("case_type", "age_months", "sex", "enrolment_season")
  missing_check <- meta_sub %>%
    select(all_of(covariate_cols)) %>%
    summarise(across(everything(), ~sum(is.na(.))))
  
  cat("\nMissing values in covariates:\n")
  print(missing_check)
  
  # Remove samples with missing covariates
  meta_sub <- meta_sub %>%
    filter(!is.na(age_months) & !is.na(sex) & !is.na(enrolment_season) & !is.na(case_type))
  
  # Update OTU table to match
  otu_sub <- otu_sub[, rownames(meta_sub)]
  
  cat("After removing samples with missing covariates:", nrow(meta_sub), "samples\n")
  
  if (nrow(meta_sub) < 10) {
    message("Skipping ", area_name, " — fewer than 10 samples after removing NAs.")
    return(NULL)
  }
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(as.matrix(otu_sub), taxa_are_rows = TRUE),
                 sample_data(meta_sub))
  
  # Set control-nonlesional_skin as the reference level
  sample_data(ps)$case_type <- factor(sample_data(ps)$case_type, 
                                      levels = c("control-nonlesional_skin", 
                                                "case-nonlesional_skin", 
                                                "case-lesional_skin"))
  
  cat("\ncase_type levels (reference first):\n")
  print(levels(sample_data(ps)$case_type))
  
  # Run ANCOM-BC2 with adjusted formula
  ancom_res <- ancombc2(
    data = ps,
    fix_formula = "case_type + age_months + sex + enrolment_season",
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "case_type",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = FALSE
  )
  
  # Extract results
  res_pairwise <- ancom_res$res_pairwise
  res_all <- ancom_res$res
  
  # Save results
  prefix <- paste0(area_name, "_ancombc2")
  
  if (is.list(res_pairwise)) {
    for (comp in names(res_pairwise)) {
      df <- res_pairwise[[comp]]
      if (is.data.frame(df)) {
        write_tsv(df, file.path(output_dir, paste0(prefix, "_pairwise_", comp, ".tsv")))
      }
    }
  }
  
  if (is.data.frame(res_all))
    write_tsv(res_all, file.path(output_dir, paste0(prefix, "_skin.tsv")))
  
  cat("Completed:", area_name, "\n")
}
###########################################################################
# Run for both Umtata and Cape Town
areas <- c("Umtata", "Cape Town")
for (loc in areas) {
  run_ancombc2_by_area(loc, metadata, otu_data, output_dir)
}
cat("\nANCOM-BC2 analysis completed for all areas.\n")