# Load required packages
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Read data
feature_table <- read.table("df_16S_filtered_feature_table_rare_Genus_absolute_10filtered_ancombc.tsv",
                            sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

metadata <- read.table("../../Metadata/differential_abundance_groups.tsv",
                       sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

metadata <- read.table(
  "../../Metadata/differential_abundance_groups.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,       # ← this sets sample IDs as row names
  check.names = FALSE
)


# Ensure sample names match
#common_samples <- intersect(rownames(feature_table), rownames(metadata))
#feature_table <- feature_table[common_samples, ]
#metadata <- metadata[common_samples, ]

# Build phyloseq object
otu <- otu_table(feature_table, taxa_are_rows = FALSE)
samp <- sample_data(metadata)
ps <- phyloseq(otu, samp)

# ✅ Run ANCOM-BC using the phyloseq object
res <- ancombc(
  data = ps,
  formula = "microbiome_type",
  group = "microbiome_type",
  p_adj_method = "fdr"
)

# Save results
write.csv(res$res$diff_abn, "diff_abundance_binary.csv")
write.csv(res$res$lfc, "log_fold_change.csv")
write.csv(res$res$p_val, "raw_pvalues.csv")
write.csv(res$res$q_val, "adjusted_qvalues.csv")
