#!/usr/bin/env python

import os
import pandas as pd
from biom import load_table, Table
from biom.util import biom_open
import logging

##########################################################################################
# SCRIPT 7: CONVERT GENUS-ASV TABLES TO RELATIVE ABUNDANCE
#           (ALL, SKIN, NASAL) FOR MULTIPLE DEPTHS
##########################################################################################

# ----------------------------------------------------------------------
# Setup logging
# ----------------------------------------------------------------------
os.makedirs('../Logs', exist_ok=True)
logging.basicConfig(
    filename='../Logs/7_relative-abundance_asv-non-collapse.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def biom_to_dataframe(biom_path):
    """Read a BIOM table and convert it to a Pandas DataFrame."""
    try:
        logging.info(f"Reading BIOM table: {biom_path}")
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe())
        logging.info(f"Loaded BIOM with shape {df.shape}")
        return df
    except Exception as e:
        logging.error(f"Error reading BIOM {biom_path}: {e}")
        return None


def convert_to_relative_abundance(df):
    """Convert a DataFrame to relative abundance (normalize each sample column)."""
    try:
        # Remove samples (columns) that sum to zero
        df_nonzero = df.loc[:, df.sum(axis=0) != 0]
        df_rel = df_nonzero.div(df_nonzero.sum(axis=0), axis=1)
        logging.info("Converted to relative abundance successfully.")
        return df_rel
    except Exception as e:
        logging.error(f"Error during relative abundance conversion: {e}")
        return None


def save_dataframe_as_biom(df, output_path):
    """Save a Pandas DataFrame as a BIOM file."""
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        biom_table = Table(df.values, observation_ids=df.index, sample_ids=df.columns)
        with biom_open(output_path, 'w') as f:
            biom_table.to_hdf5(f, generated_by="Relative Abundance Script")
        logging.info(f"Saved BIOM: {output_path}")
    except Exception as e:
        logging.error(f"Error saving BIOM {output_path}: {e}")

# ----------------------------------------------------------------------
# Main Execution
# ----------------------------------------------------------------------
if __name__ == "__main__":
    prevalence_thresholds = ['10pct', '5pct', '1pct', '0pct']
    rarefaction_depths = [350, 1000, 1500, 2000]  # Based on QIIME2 rarefaction curve
    specimens = ['all', 'skin', 'nasal']

    input_dir = "../Data/Tables/Count_Tables"
    output_dir = "../Data/Tables/Relative_Abundance_Tables"
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Starting relative abundance generation for Genus-ASV tables (ALL, SKIN, NASAL).")

    for threshold in prevalence_thresholds:
        for depth in rarefaction_depths:
            for specimen in specimens:
                biom_path = os.path.join(
                    input_dir,
                    f"6_209766_feature_table_dedup_prev-filt-{threshold}_rare-{depth}_Genus-ASV_{specimen}.biom"
                )

                if not os.path.exists(biom_path):
                    logging.warning(f"BIOM not found: {biom_path}")
                    continue

                print(f"\n=== Processing prevalence {threshold}, depth {depth}, specimen {specimen} ===")

                # Step 1: Load BIOM
                df = biom_to_dataframe(biom_path)
                if df is None or df.empty:
                    logging.warning(f"Skipping {biom_path} — DataFrame empty or failed to load.")
                    continue

                # Step 2: Convert to relative abundance
                df_rel = convert_to_relative_abundance(df)
                if df_rel is None or df_rel.empty:
                    logging.warning(f"Relative abundance conversion failed for {specimen}")
                    continue

                # Step 3: Save output BIOM
                output_file = os.path.join(
                    output_dir,
                    f"7_209766_feature_table_dedup_prev-filt-{threshold}_rare-{depth}_Genus-ASV_{specimen}_rel.biom"
                )
                save_dataframe_as_biom(df_rel, output_file)
                print(f"  → Saved relative abundance BIOM: {output_file}")

    logging.info("Completed relative abundance table creation for all specimen subsets.")
    print("Done. Log written to: ../Logs/7_relative-abundance_asv-non-collapse.log")
