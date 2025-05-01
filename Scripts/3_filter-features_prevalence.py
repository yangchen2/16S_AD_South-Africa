#!/usr/bin/env python

import pandas as pd
from biom import load_table, Table
from biom.util import biom_open
import os
import logging
import matplotlib.pyplot as plt

# Ensure required directories exist
os.makedirs('../Logs', exist_ok=True)
os.makedirs('../Plots/Dataset_stats', exist_ok=True)

# Setup logging
logging.basicConfig(
    filename='../Logs/3_filter_features_prevalence.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def read_biom_to_df(biom_path):
    logging.info(f"Loading BIOM file: {biom_path}")
    table = load_table(biom_path)
    df = pd.DataFrame(table.to_dataframe(dense=True)).T  # Samples as rows
    logging.info(f"BIOM loaded with shape: {df.shape}")
    return df

def calculate_prevalence(df):
    n_samples = df.shape[0]
    prevalence = (df > 0).sum(axis=0) / n_samples * 100  # convert to percentage
    logging.info(f"Calculated prevalence for {len(prevalence)} features.")
    return prevalence

def plot_prevalence(prevalence_series, output_path):
    plt.figure(figsize=(8, 10))
    plt.hist(prevalence_series, bins=50, color='steelblue', edgecolor='black')
    plt.title('Feature Prevalence Distribution')
    plt.xlabel('Prevalence (% of samples)')
    plt.ylabel('Number of Features')
    plt.gca().invert_yaxis()  # High prevalence at top
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    logging.info(f"Prevalence plot saved to: {output_path}")

def filter_by_prevalence(df, prevalence_series, threshold):
    filtered_features = prevalence_series[prevalence_series >= threshold].index
    filtered_df = df[filtered_features]
    logging.info(f"Features before filtering: {df.shape[1]}")
    logging.info(f"Features after filtering at {threshold}%: {filtered_df.shape[1]}")
    return filtered_df

def save_df_as_biom(df, output_path):
    logging.info(f"Saving filtered table to: {output_path}")
    df_T = df.T  # BIOM format expects features as rows
    biom_table = Table(df_T.values, observation_ids=df_T.index, sample_ids=df_T.columns)
    with biom_open(output_path, 'w') as f:
        biom_table.to_hdf5(f, "Filtered by prevalence")
    logging.info(f"Filtered BIOM saved successfully.")

if __name__ == '__main__':
    try:
        # HARD-CODED input biom
        biom_path = '../Data/Tables/Absolute_Abundance_Tables/209766_filtered_feature_table.biom'

        # Load the raw data
        df = read_biom_to_df(biom_path)
        prevalence = calculate_prevalence(df)

        # Plot raw prevalence distribution once
        plot_prevalence(prevalence, '../Plots/Dataset_stats/feature_prevalence_distribution.png')

        # Loop through thresholds
        thresholds = [10, 5, 1, 0]  # percent thresholds
        for threshold in thresholds:
            output_path = f'../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_{threshold}pct.biom'
            filtered_df = filter_by_prevalence(df, prevalence, threshold=threshold)
            save_df_as_biom(filtered_df, output_path)
            logging.info(f"Completed filtering at {threshold}% threshold.")

        logging.info("Filtering and plotting completed successfully for all thresholds.")
        print("Done. Log written to: ../Logs/3_filter_features_prevalence.log")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
