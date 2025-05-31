#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import warnings
import logging

warnings.simplefilter(action='ignore', category=FutureWarning)

# Setup logging
os.makedirs('../logs', exist_ok=True)
logging.basicConfig(
    filename='../logs/6_taxonomy_tbl_map_ASV-non-collapse.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Read in taxonomy
logging.info("Loading GreenGreens2 taxonomy...")
gg_taxonomy = q2.Artifact.load('../Reference/2022.10.taxonomy.asv.tsv.qza').view(pd.DataFrame)
logging.info("Taxonomy loaded")

def add_unique_tax_labels(tbl_path: str, level: str, prevalence: str):
    """
    Map ASV IDs to taxonomy names (e.g., g__Streptococcus_ASV-1) ordered by total abundance,
    without collapsing rows.
    Output one biom file per taxonomy level and prevalence threshold.
    """
    # Load the BIOM table
    biom_table = load_table(tbl_path)
    df = pd.DataFrame(biom_table.to_dataframe()).transpose()

    # Attach taxonomy
    df_tax = df.transpose()
    df_tax = df_tax.merge(gg_taxonomy['Taxon'], how='left', left_index=True, right_index=True)

    if df_tax['Taxon'].isnull().any():
        logging.warning(f"Some ASVs are missing taxonomy assignments at {level} level for {prevalence}.")

    # Split taxonomy into levels
    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_tax[taxonomy_levels] = df_tax['Taxon'].str.split(';', expand=True)

    # Drop ASVs without taxonomy at the desired level
    df_tax = df_tax[df_tax[level].notnull()]

    # Compute total abundance across all samples
    df_tax['TotalAbundance'] = df_tax.drop(columns=taxonomy_levels + ['Taxon']).sum(axis=1)

    # Assign new names ordered by abundance within each taxon
    new_index = []
    for taxon, group in df_tax.groupby(level):
        group_sorted = group.sort_values('TotalAbundance', ascending=False)
        for i, idx in enumerate(group_sorted.index, 1):
            clean_taxon = taxon.strip().replace(' ', '_')  # Clean for filenames
            new_id = f"{clean_taxon}_ASV-{i}"
            new_index.append((idx, new_id))

    # Apply new index
    index_map = dict(new_index)
    df_tax.rename(index=index_map, inplace=True)

    # Drop helper columns
    df_tax.drop(columns=['TotalAbundance'], inplace=True)
    df_out = df_tax.drop(columns=taxonomy_levels + ['Taxon'])

    # Save new BIOM table
    os.makedirs("../Data/Tables/Absolute_Abundance_Tables/", exist_ok=True)
    obs_ids = df_out.index
    samp_ids = df_out.columns
    biom_out = biom.table.Table(df_out.values, observation_ids=obs_ids, sample_ids=samp_ids)

    # out_file = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare_{level}-ASV-non-collapse.biom"
    out_file = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare_filtered_{level}-ASV-non-collapse.biom" # for alpha diversity

    with biom_open(out_file, 'w') as f:
        biom_out.to_hdf5(f, generated_by=f"Taxonomy mapping at {level} level, {prevalence}")

    logging.info(f"Saved mapped BIOM file for {level} at {prevalence}: {out_file}")

if __name__ == '__main__':
    try:
        prevalence_thresholds = ['1pct']
        taxonomy_level = "Genus"

        for prevalence in prevalence_thresholds:
            # biom_path = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare.biom"
            biom_path = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare_filtered.biom" # for alpha diversity
            
            if not os.path.exists(biom_path):
                logging.warning(f"BIOM file not found for {prevalence} threshold: {biom_path}")
                continue

            logging.info(f"Processing {prevalence} threshold at {taxonomy_level} level")
            add_unique_tax_labels(biom_path, level=taxonomy_level, prevalence=prevalence)

        logging.info("All taxonomy mappings completed for all prevalence thresholds.")
        print("Done. Log written to: ../logs/7_taxonomy_tbl_map_ASV-non-collapse.log")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise