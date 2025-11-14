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
#           RAREFIED: collapse to Genus-ASV-X
#           NON-RAREFIED: ASV-non-collapsed (original ASVs)
##########################################################################################

warnings.simplefilter(action='ignore', category=FutureWarning)
WRITE_FASTA = True
WRITE_BIOM = True

# ------------------------------------------------------------------
# Logging setup
# ------------------------------------------------------------------
os.makedirs("../logs", exist_ok=True)
logging.basicConfig(
    filename="../logs/6_taxonomy_tbl_map_Genus_ASV_with-fasta.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ------------------------------------------------------------------
# Load Greengenes2 taxonomy
# ------------------------------------------------------------------
logging.info("Loading Greengenes2 taxonomy artifact...")
gg_taxonomy = q2.Artifact.load("../Reference/2022.10.taxonomy.asv.tsv.qza").view(pd.DataFrame)
logging.info("Greengenes2 taxonomy successfully loaded.")

# ------------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------------
metadata_path = "../Metadata/16S_AD_South-Africa_metadata_subset.tsv"
metadata = pd.read_csv(metadata_path, sep="\t").set_index("#sample-id")

# ------------------------------------------------------------------
# FASTA writer
# ------------------------------------------------------------------
def write_fasta_from_index_map(index_map: dict, output_fasta_path: str):
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


# ------------------------------------------------------------------
# Genus-ASV collapsing for rarefied tables
# ------------------------------------------------------------------
def add_unique_tax_labels(tbl_path: str, level: str):
    biom_table = load_table(tbl_path)
    df = pd.DataFrame(biom_table.to_dataframe()).T  # samples × ASVs

    df_tax = df.T.merge(gg_taxonomy["Taxon"], how="left", left_index=True, right_index=True)
    taxonomy_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    df_tax["Taxon"] = df_tax["Taxon"].fillna("Unknown")
    df_tax[taxonomy_levels] = df_tax["Taxon"].str.split(";", expand=True)
    df_tax = df_tax[df_tax[level].notnull()]

    df_tax["TotalAbundance"] = df_tax.drop(columns=taxonomy_levels + ["Taxon"]).sum(axis=1)

    index_map = {}
    for taxon, group in df_tax.groupby(level):
        sorted_group = group.sort_values("TotalAbundance", ascending=False)
        clean_taxon = taxon.strip().replace(" ", "_").replace(":", "_")
        for i, idx in enumerate(sorted_group.index, 1):
            index_map[idx] = f"{clean_taxon}_ASV-{i}"

    df_tax_renamed = df_tax.rename(index=index_map)
    df_out = df_tax_renamed.drop(columns=taxonomy_levels + ["Taxon", "TotalAbundance"])

    return index_map, df_out.T


# ------------------------------------------------------------------
# Specimen filtering
# ------------------------------------------------------------------
def filter_samples_by_specimen(df_counts: pd.DataFrame, specimen: str):
    if specimen is None:
        return df_counts
    valid_ids = metadata[metadata["specimen"] == specimen].index
    return df_counts.loc[df_counts.index.intersection(valid_ids)]


# ------------------------------------------------------------------
# BIOM saver
# ------------------------------------------------------------------
def save_biom_table(df_counts: pd.DataFrame, output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    biom_out = biom.table.Table(df_counts.T.values, observation_ids=df_counts.columns, sample_ids=df_counts.index)
    with biom_open(output_path, "w") as f:
        biom_out.to_hdf5(f, generated_by="ASV + Genus-ASV mapping")
    logging.info(f"Saved BIOM: {output_path}")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
if __name__ == "__main__":
    try:
        prevalence_thresholds = ["10pct", "5pct", "1pct", "0pct"]
        taxonomy_level = "Genus"
        rare_depths = [350, 1000, 1500, 2000]

        # rarefied (collapse)
        biom_prefix_rarefied = ["5_"]

        # non-rarefied (ASV-non-collapsed)
        biom_prefix_nonraref = ["3_"]

        variants = {"all": None, "skin": "skin", "nasal": "nasal"}

        # =============================================================
        # 1. RAREFIED: Genus-ASV collapsed
        # =============================================================
        for prevalence in prevalence_thresholds:
            for depth in rare_depths:
                for biom_prefix in biom_prefix_rarefied:
                    biom_path = (
                        f"../Data/Tables/Count_Tables/"
                        f"{biom_prefix}209766_feature_table_dedup_prev-filt-{prevalence}_rare-{depth}.biom"
                    )

                    if not os.path.exists(biom_path):
                        continue

                    print(f"\n=== RAREFIED: prevalence={prevalence}, depth={depth} ===")

                    index_map, df_counts = add_unique_tax_labels(biom_path, taxonomy_level)

                    for variant, specimen in variants.items():
                        df_sub = filter_samples_by_specimen(df_counts, specimen)
                        df_sub = df_sub.loc[:, df_sub.sum(axis=0) > 0]

                        # BIOM output
                        biom_out = (
                            f"../Data/Tables/Count_Tables/"
                            f"6_209766_feature_table_dedup_prev-filt-{prevalence}_rare-{depth}_Genus-ASV_{variant}.biom"
                        )
                        save_biom_table(df_sub, biom_out)

                        # FASTA output
                        fasta_out = (
                            f"../Data/Fasta/"
                            f"rare-{depth}_prev-{prevalence}_Genus-ASV_{variant}.fasta"
                        )
                        index_sub = {k: v for k, v in index_map.items() if v in df_sub.columns}
                        write_fasta_from_index_map(index_sub, fasta_out)

        # =============================================================
        # 2. NON-RAREFIED: ASV-non-collapsed (original ASVs)
        # =============================================================
        for prevalence in prevalence_thresholds:
            for biom_prefix in biom_prefix_nonraref:
                biom_path_nr = (
                    f"../Data/Tables/Count_Tables/"
                    f"{biom_prefix}209766_feature_table_dedup_prev-filt-{prevalence}.biom"
                )

                if not os.path.exists(biom_path_nr):
                    continue

                print(f"\n=== NON-RAREFIED: ASV-non-collapsed, prevalence={prevalence} ===")

                biom_table_nr = load_table(biom_path_nr)
                df_nr = pd.DataFrame(biom_table_nr.to_dataframe()).T

                for variant, specimen in variants.items():
                    df_nr_sub = filter_samples_by_specimen(df_nr, specimen)
                    df_nr_sub = df_nr_sub.loc[:, df_nr_sub.sum(axis=0) > 0]

                    # BIOM output
                    biom_out_nr = (
                        f"../Data/Tables/Count_Tables/"
                        f"6_209766_feature_table_dedup_prev-filt-{prevalence}_Genus-ASV_{variant}.biom"
                    )
                    save_biom_table(df_nr_sub, biom_out_nr)


        logging.info("All processing complete.")
        print("Done. Log written to ../logs/6_taxonomy_tbl_map_Genus_ASV_with-fasta.log")

    except Exception as e:
        logging.error(f"Error occurred: {e}")
        raise
