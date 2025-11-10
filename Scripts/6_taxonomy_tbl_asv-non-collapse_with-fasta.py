#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import warnings
import logging
from skbio import DNA
from skbio.io import write
import qiime2 as q2

##########################################################################################
# SCRIPT 6: MAP ASV SEQUENCES TO TAXONOMY NAMES AND EXPORT FASTA + BIOM
#           (ALL, SKIN, NASAL) FOR MULTIPLE DEPTHS — ENSURES FEATURE MATCHING
##########################################################################################

warnings.simplefilter(action='ignore', category=FutureWarning)
WRITE_FASTA = True
WRITE_BIOM = True

# ------------------------------------------------------------------
# Setup logging
# ------------------------------------------------------------------
os.makedirs("../logs", exist_ok=True)
logging.basicConfig(
    filename="../logs/6_taxonomy_tbl_map_Genus_ASV_with-fasta.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ------------------------------------------------------------------
# Load taxonomy + metadata
# ------------------------------------------------------------------
logging.info("Loading Greengenes2 taxonomy artifact...")
gg_taxonomy = q2.Artifact.load("../Reference/2022.10.taxonomy.asv.tsv.qza").view(pd.DataFrame)
logging.info("Greengenes2 taxonomy successfully loaded.")

metadata_path = "../Metadata/16S_AD_South-Africa_metadata_subset.tsv"
metadata = pd.read_csv(metadata_path, sep="\t")
metadata.set_index("#sample-id", inplace=True)

# ------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------
def write_fasta_from_index_map(index_map: dict, output_fasta_path: str):
    """Write FASTA where headers are taxonomy-labeled ASV IDs."""
    os.makedirs(os.path.dirname(output_fasta_path), exist_ok=True)
    count = 0
    with open(output_fasta_path, "w") as f:
        for seq, label in index_map.items():
            try:
                record = DNA(seq.strip(), metadata={"id": label})
                write(record, format="fasta", into=f)
                count += 1
            except Exception as e:
                logging.warning(f"Skipped invalid sequence: {seq[:15]}... ({e})")
    logging.info(f"FASTA written: {output_fasta_path} ({count} sequences)")
    print(f"  → {count} sequences written to {output_fasta_path}")


def add_unique_tax_labels(tbl_path: str, level: str, prevalence: str, depth: int):
    """Map ASV IDs (ASV sequences) to taxonomy-labeled IDs like g__Streptococcus_ASV-1."""
    biom_table = load_table(tbl_path)
    df = pd.DataFrame(biom_table.to_dataframe()).T  # samples × features
    logging.info(f"Loaded BIOM {tbl_path} with shape {df.shape}")

    # Attach taxonomy
    df_tax = df.T
    df_tax = df_tax.merge(gg_taxonomy["Taxon"], how="left", left_index=True, right_index=True)
    df_tax["Taxon"] = df_tax["Taxon"].fillna("Unknown")

    taxonomy_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    df_tax[taxonomy_levels] = df_tax["Taxon"].str.split(";", expand=True)
    df_tax = df_tax[df_tax[level].notnull()]

    # Compute total abundance & assign ASV labels
    df_tax["TotalAbundance"] = df_tax.drop(columns=taxonomy_levels + ["Taxon"]).sum(axis=1)
    index_map = {}
    for taxon, group in df_tax.groupby(level):
        sorted_group = group.sort_values("TotalAbundance", ascending=False)
        clean_taxon = taxon.strip().replace(" ", "_").replace(":", "_")
        for i, idx in enumerate(sorted_group.index, 1):
            index_map[idx] = f"{clean_taxon}_ASV-{i}"

    # Replace feature IDs
    df_tax_renamed = df_tax.rename(index=index_map)
    df_out = df_tax_renamed.drop(columns=taxonomy_levels + ["Taxon", "TotalAbundance"])

    return index_map, df_out.T  # sample × feature table


def filter_samples_by_specimen(df_counts: pd.DataFrame, specimen: str) -> pd.DataFrame:
    """Subset BIOM DataFrame to samples of a given specimen type."""
    if specimen is None:
        return df_counts
    valid_samples = metadata[metadata["specimen"] == specimen].index
    subset = df_counts.loc[df_counts.index.intersection(valid_samples)]
    if subset.empty:
        logging.warning(f"No samples found for specimen={specimen}")
    return subset


def save_biom_table(df_counts: pd.DataFrame, output_path: str):
    """Save BIOM table."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    biom_out = biom.table.Table(df_counts.T.values, observation_ids=df_counts.columns, sample_ids=df_counts.index)
    with biom_open(output_path, "w") as f:
        biom_out.to_hdf5(f, generated_by="ASV mapping per specimen type")
    logging.info(f"Saved BIOM: {output_path}")

# ------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------
if __name__ == "__main__":
    try:
        prevalence_thresholds = ["10pct", "5pct", "1pct", "0pct"]
        taxonomy_level = "Genus"
        rarefaction_depths = [350, 1000, 1500, 2000]
        biom_prefixes = ["5_", "3_"]

        for prevalence in prevalence_thresholds:
            for depth in rarefaction_depths:
                for biom_prefix in biom_prefixes:
                    biom_path = (
                        f"../Data/Tables/Count_Tables/"
                        f"{biom_prefix}209766_feature_table_dedup_prev-filt-{prevalence}_rare-{depth}.biom"
                    )

                    if not os.path.exists(biom_path):
                        logging.warning(f"BIOM not found: {biom_path}")
                        continue

                    print(f"\n=== Processing prefix {biom_prefix}, prevalence {prevalence}, depth {depth} ===")
                    index_map, df_counts_all = add_unique_tax_labels(biom_path, taxonomy_level, prevalence, depth)

                    variants = {
                        "all": None,
                        "skin": "skin",
                        "nasal": "nasal"
                    }

                    for variant, specimen in variants.items():
                        df_counts_subset = filter_samples_by_specimen(df_counts_all, specimen)
                        if df_counts_subset.empty:
                            print(f"No samples found for {variant} at depth {depth}. Skipping.")
                            continue

                        features_present = df_counts_subset.columns[df_counts_subset.sum(axis=0) > 0].tolist()
                        df_counts_subset = df_counts_subset[features_present]
                        index_map_subset = {seq: label for seq, label in index_map.items() if label in features_present}

                        n_biom = len(df_counts_subset.columns)
                        n_fasta = len(index_map_subset)
                        if n_biom != n_fasta:
                            logging.warning(f"Feature mismatch for {variant} (BIOM={n_biom}, FASTA={n_fasta}) — syncing.")
                            shared_features = set(features_present).intersection(index_map_subset.values())
                            df_counts_subset = df_counts_subset[list(shared_features)]
                            index_map_subset = {k: v for k, v in index_map_subset.items() if v in shared_features}
                            n_biom = len(df_counts_subset.columns)
                            n_fasta = len(index_map_subset)

                        logging.info(f"Saving {variant}: {n_biom} features, depth={depth}")

                        # Save BIOM
                        if WRITE_BIOM:
                            biom_outfile = (
                                f"../Data/Tables/Count_Tables/"
                                f"6_{biom_prefix}209766_feature_table_dedup_prev-filt-{prevalence}_rare-{depth}_Genus-ASV_{variant}.biom"
                            )
                            save_biom_table(df_counts_subset, biom_outfile)

                        # Save FASTA
                        if WRITE_FASTA:
                            fasta_out = (
                                f"../Data/Fasta/"
                                f"{biom_prefix}209766_feature_table_dedup_prev-filt-{prevalence}_rare-{depth}_Genus-ASV_{variant}.fasta"
                            )
                            write_fasta_from_index_map(index_map_subset, fasta_out)

                        assert n_biom == n_fasta, (
                            f"Mismatch detected after filtering for {variant} "
                            f"(BIOM={n_biom}, FASTA={n_fasta})"
                        )

                    logging.info(f"Completed: prefix={biom_prefix}, prevalence={prevalence}, depth={depth}")

        logging.info("All taxonomy mappings, BIOMs, and FASTAs completed (ALL, SKIN, NASAL).")
        print("Done. Log written to: ../logs/6_taxonomy_tbl_map_Genus_ASV_with-fasta.log")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise
