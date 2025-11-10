#!/usr/bin/env python

import os
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
from numpy.random import RandomState
import matplotlib.pyplot as plt
import logging

##########################################################################################
# SCRIPT 5: PERFORMS RAREFACTION ON FEATURE TABLE
##########################################################################################

# Setup logging
os.makedirs('../Logs', exist_ok=True)
logging.basicConfig(filename='../Logs/5_rarefaction.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def read_and_convert_biom(biom_path: str) -> pd.DataFrame:
    logging.info(f"Loading BIOM file: {biom_path}")
    try:
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe()).T
        df.index = df.index.str.replace('^15564\.', '', regex=True)
        logging.info(f"BIOM table shape: {df.shape}")
        return df
    except Exception as e:
        logging.error(f"Error in processing BIOM file: {e}")
        raise


def load_metadata(metadata_path: str) -> pd.DataFrame:
    try:
        metadata = pd.read_csv(metadata_path, sep='\t')
        metadata.set_index('#sample-id', inplace=True)
        logging.info(f"Metadata loaded with shape: {metadata.shape}")
        return metadata
    except Exception as e:
        logging.error(f"Error in loading metadata: {e}")
        raise


def rarefy_table(biom_df: pd.DataFrame, seed: int = 42, depth: int = None) -> pd.DataFrame:
    if depth is None:
        depth = int(biom_df.sum(axis=1).min())
        logging.info(f"Min depth: {depth}")

    if depth is None:
        depth = int(biom_df.sum(axis=1).min())
        logging.info(f"Mean depth: {depth}")

    logging.info(f"Using rarefaction depth: {depth}")

    prng = RandomState(seed)
    rarefied_data = np.zeros_like(biom_df.values)

    for i, row in enumerate(biom_df.values):
        total_count = row.sum()
        if total_count >= depth:
            probabilities = row / total_count
            chosen_indices = prng.choice(biom_df.columns.size, depth, p=probabilities)
            rarefied_row = np.bincount(chosen_indices, minlength=biom_df.columns.size)
            rarefied_data[i] = rarefied_row

    rarefied_df = pd.DataFrame(rarefied_data, index=biom_df.index, columns=biom_df.columns)
    logging.info("Rarefaction completed.")
    return rarefied_df


def save_as_biom(df: pd.DataFrame, output_path: str):
    df = df.transpose()
    table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with biom_open(output_path, 'w') as f:
        table.to_hdf5(f, "rarefaction script")
    logging.info(f"Saved rarefied BIOM table to {output_path}")



if __name__ == '__main__':
    try:
        prevalence_thresholds = ['10pct', '5pct', '1pct', '0pct']
        metadata_path = '../Metadata/16S_AD_South-Africa_metadata_subset.tsv'
        
        chosen_depths = [350, 1000, 1500, 2000]  # Based on QIIME2 rarefaction curve
        
        for chosen_depth in chosen_depths:
            logging.info(f"Chosen rarefaction depth: {chosen_depth}")
            
            for prevalence in prevalence_thresholds:
                biom_path = f'../Data/Tables/Count_Tables/3_209766_feature_table_dedup_prev-filt-{prevalence}.biom'
                filtered_output_path = f'../Data/Tables/Count_Tables/5_209766_feature_table_dedup_prev-filt-{prevalence}_rare-{chosen_depth}.biom'

                if not os.path.exists(biom_path):
                    logging.warning(f"BIOM file not found for {prevalence}: {biom_path}")
                    continue

                # Load BIOM + metadata
                df = read_and_convert_biom(biom_path)
                metadata = load_metadata(metadata_path)

                # Filter out samples below that rarefaction depth
                df_filtered = df[df.sum(axis=1) >= chosen_depth]
                logging.info(f"Samples retained after filtering: {df_filtered.shape[0]}")

                # Perform rarefaction
                rarefied_filtered_df = rarefy_table(df_filtered, seed=42, depth=chosen_depth)
                save_as_biom(rarefied_filtered_df, filtered_output_path)
                logging.info(f"Saved rarefied BIOM table with filtered samples to {filtered_output_path}")

            logging.info(f"Rarefaction completed for depth {chosen_depth}")

        logging.info("Rarefaction pipeline completed for all prevalence thresholds and depths.")
        print("Done. Log written to: ../Logs/5_rarefaction.log")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

