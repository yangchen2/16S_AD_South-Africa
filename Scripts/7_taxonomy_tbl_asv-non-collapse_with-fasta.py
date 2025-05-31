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

warnings.simplefilter(action='ignore', category=FutureWarning)

# Toggle FASTA output
WRITE_FASTA = True

# Setup logging
os.makedirs('../logs', exist_ok=True)
logging.basicConfig(
    filename='../logs/7_taxonomy_tbl_map_ASV-non-collapse.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Read in taxonomy
logging.info("Loading GreenGreens2 taxonomy...")
import qiime2 as q2
gg_taxonomy = q2.Artifact.load('../Reference/2022.10.taxonomy.asv.tsv.qza').view(pd.DataFrame)
logging.info("Taxonomy loaded")

def write_fasta_with_taxonomy_from_columns(index_map: dict, output_fasta_path: str):
    """
    Write a FASTA file where each header is the taxonomy-labeled ASV name (e.g., g__Streptococcus_ASV-1)
    and the sequence is taken directly from the ASV ID (i.e., column name).
    """
    logging.info("Writing FASTA using ASV sequences as column names")

    with open(output_fasta_path, 'w') as f:
        for old_id, new_id in index_map.items():
            sequence = old_id  # ASV sequence is the column name
            record = DNA(sequence, metadata={'id': new_id})
            write(record, format='fasta', into=f)

    logging.info(f"FASTA file written to {output_fasta_path}")

def add_unique_tax_labels(tbl_path: str, level: str, prevalence: str):
    """
    Map ASV IDs (ASV sequences as column names) to taxonomy names (e.g., g__Streptococcus_ASV-1),
    ordered by total abundance, without collapsing rows.
    Returns:
        index_map: dict mapping original ASV IDs to taxonomy-labeled IDs
    """
    # Load the BIOM table
    biom_table = load_table(tbl_path)
    df = pd.DataFrame(biom_table.to_dataframe()).transpose()
    print('Table shape: ' + str(df.shape))

    # Attach taxonomy
    df_tax = df.transpose()
    df_tax = df_tax.merge(gg_taxonomy['Taxon'], how='left', left_index=True, right_index=True)

    if df_tax['Taxon'].isnull().any():
        #logging.warning(f"Some ASVs are missing taxonomy assignments at {level} level for {prevalence}.")
        df_tax['Taxon'] = df_tax['Taxon'].fillna('Unknown')

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
            clean_taxon = taxon.strip().replace(' ', '_')
            new_id = f"{clean_taxon}_ASV-{i}"
            new_index.append((idx, new_id))

    index_map = dict(new_index)
    df_tax_renamed = df_tax.rename(index=index_map)

    # Drop helper columns
    df_tax_renamed.drop(columns=['TotalAbundance'], inplace=True)
    df_out = df_tax_renamed.drop(columns=taxonomy_levels + ['Taxon'])

    # Save new BIOM table
    os.makedirs("../Data/Tables/Absolute_Abundance_Tables/", exist_ok=True)
    obs_ids = df_out.index
    samp_ids = df_out.columns
    biom_out = biom.table.Table(df_out.values, observation_ids=obs_ids, sample_ids=samp_ids)

    out_file = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare_{level}-ASV-non-collapse.biom"
    with biom_open(out_file, 'w') as f:
        biom_out.to_hdf5(f, generated_by=f"Taxonomy mapping at {level} level, {prevalence}")

    logging.info(f"Saved mapped BIOM file for {level} at {prevalence}: {out_file}")
    return index_map

if __name__ == '__main__':
    try:
        prevalence_thresholds = ['10pct', '5pct', '1pct', '0pct']
        taxonomy_level = "Genus"

        for prevalence in prevalence_thresholds:
            biom_path = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}.biom"
            if not os.path.exists(biom_path):
                logging.warning(f"BIOM file not found for {prevalence} threshold: {biom_path}")
                continue

            logging.info(f"Processing {prevalence} threshold at {taxonomy_level} level")
            index_map = add_unique_tax_labels(biom_path, level=taxonomy_level, prevalence=prevalence)

            if WRITE_FASTA:
                fasta_out_path = f"../Data/Fasta/209766_filtered_by_prevalence_{prevalence}.fasta"
                write_fasta_with_taxonomy_from_columns(index_map, fasta_out_path)

        logging.info("All taxonomy mappings and FASTA exports completed.")
        print("Done. Log written to: ../logs/7_taxonomy_tbl_map_ASV-non-collapse.log")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise
