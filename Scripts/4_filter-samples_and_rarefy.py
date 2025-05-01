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

# Setup logging
os.makedirs('../Logs', exist_ok=True)
logging.basicConfig(filename='../Logs/4_filter-samples_and_rarefy.log', level=logging.INFO,
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

def plot_sorted_read_counts(biom_df: pd.DataFrame, metadata: pd.DataFrame, prevalence: str):
    sample_read_counts = biom_df.sum(axis=1).sort_values(ascending=False)
    sample_read_counts_df = sample_read_counts.to_frame(name='read_count').merge(
        metadata[['case_type']], left_index=True, right_index=True, how='left'
    ).dropna(subset=['case_type']).reset_index(drop=True)

    color_map = {
        'control-nonlesional skin': 'blue',
        'case-nonlesional skin': 'orange',
        'case-lesional skin': 'red',
        'control-anterior nares': 'green',
        'case-anterior nares': 'pink'
    }
    sample_read_counts_df['color'] = sample_read_counts_df['case_type'].map(color_map)

    plt.figure(figsize=(10, 6))
    plt.scatter(x=sample_read_counts_df.index, y=sample_read_counts_df['read_count'],
                c=sample_read_counts_df['color'], s=20)

    legend_order = list(color_map.keys())
    for case_type in legend_order:
        plt.scatter([], [], color=color_map[case_type], label=case_type)

    plt.title(f'Read Counts per Sample (Sorted High to Low) - {prevalence}')
    plt.xlabel('Sample (sorted by read count)')
    plt.ylabel('Read Count')
    plt.xticks([])
    plt.legend(title="Case Type")
    plt.grid(True)
    os.makedirs('../Plots/Dataset_stats/', exist_ok=True)
    plt.savefig(f'../Plots/Dataset_stats/read_counts_by_case_type_{prevalence}.png', dpi=600)
    logging.info(f"Read counts plot saved for {prevalence}.")

def plot_rarefaction_depths(biom_df: pd.DataFrame, prevalence: str) -> int:
    sample_depths = biom_df.sum(axis=1)
    mean_depth = int(sample_depths.mean())

    plt.figure(figsize=(10, 6))
    plt.hist(sample_depths, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(sample_depths.min(), color='red', linestyle='--', label=f'Min Depth: {sample_depths.min()}')
    plt.axvline(mean_depth, color='blue', linestyle='--', label=f'Mean Depth: {mean_depth}')
    plt.title(f'Distribution of Sample Read Depths - {prevalence}')
    plt.xlabel('Read Depth')
    plt.ylabel('Frequency')
    plt.legend()
    os.makedirs('../Plots/Dataset_stats/', exist_ok=True)
    plt.savefig(f'../Plots/Dataset_stats/read_counts_by_sample_{prevalence}.png', dpi=600)
    logging.info(f"Rarefaction depth plot saved for {prevalence}. Mean depth: {mean_depth}")
    return mean_depth

def rarefy_table(biom_df: pd.DataFrame, seed: int = 42, depth: int = None) -> pd.DataFrame:
    if depth is None:
        depth = int(biom_df.sum(axis=1).min())
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
        metadata_path = '../Data/Metadata/updated_clean_ant_skin_metadata.tab'

        for prevalence in prevalence_thresholds:
            biom_path = f'../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}.biom'
            output_path = f'../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{prevalence}_rare.biom'

            if not os.path.exists(biom_path):
                logging.warning(f"BIOM file not found for {prevalence}: {biom_path}")
                continue

            df = read_and_convert_biom(biom_path)
            metadata = load_metadata(metadata_path)

            plot_sorted_read_counts(df, metadata, prevalence)
            mean_depth = plot_rarefaction_depths(df, prevalence)

            # You can adjust the rarefaction depth if you want based on your project:
            rarefied_df = rarefy_table(df, seed=42, depth=350)

            save_as_biom(rarefied_df, output_path)

        logging.info("Rarefaction pipeline completed for all prevalence thresholds.")
        print("Done. Log written to: ../Logs/4_filter-samples_and_rarefy.log")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
