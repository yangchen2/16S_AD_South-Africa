#!/usr/bin/env python

import pandas as pd
from biom import load_table, Table
from biom.util import biom_open
import qiime2 as q2
import os
import logging
import matplotlib.pyplot as plt

##########################################################################################
# SCRIPT 3: FILTERS FEATURES BASED ON SAMPLE PREVALENCE
##########################################################################################

# Ensure required directories exist
os.makedirs('../Logs', exist_ok=True)
os.makedirs('../Data/Tables/Count_Tables', exist_ok=True)

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


def filter_by_prevalence(df, prevalence_series, threshold):
    filtered_features = prevalence_series[prevalence_series >= threshold].index
    filtered_df = df[filtered_features]
    logging.info(f"Features before filtering: {df.shape[1]}")
    logging.info(f"Features after filtering at {threshold}%: {filtered_df.shape[1]}")
    return filtered_df

def save_df_as_biom_and_qza(df, biom_path):
    """
    Save filtered DataFrame as both BIOM and QIIME2 Artifact (.qza)
    """
    logging.info(f"Saving filtered table to: {biom_path}")
    df_T = df.T  # BIOM expects features as rows
    biom_table = Table(df_T.values, observation_ids=df_T.index, sample_ids=df_T.columns)
    with biom_open(biom_path, 'w') as f:
        biom_table.to_hdf5(f, "Filtered by prevalence")
    logging.info("Filtered BIOM saved successfully.")

    # Convert BIOM to QZA
    try:
        qza_path = biom_path.replace('.biom', '.qza')
        artifact = q2.Artifact.import_data('FeatureTable[Frequency]', biom_path)
        artifact.save(qza_path)
        logging.info(f"QIIME2 artifact saved: {qza_path}")
    except Exception as e:
        logging.error(f"Error converting to QZA: {e}")
        raise

if __name__ == '__main__':
    try:
        biom_path = '../Data/Tables/Count_Tables/2_209766_feature_table_dedup.biom'

        # Load BIOM and compute prevalence
        df = read_biom_to_df(biom_path)
        prevalence = calculate_prevalence(df)

        # Loop through prevalence thresholds
        thresholds = [10, 5, 1, 0]
        for threshold in thresholds:
            biom_out = f'../Data/Tables/Count_Tables/3_209766_feature_table_dedup_prev-filt-{threshold}pct.biom'
            filtered_df = filter_by_prevalence(df, prevalence, threshold)
            save_df_as_biom_and_qza(filtered_df, biom_out)
            logging.info(f"Completed filtering at {threshold}% threshold.")

        logging.info("Filtering and QZA export completed successfully for all thresholds.")
        print("Done. Log written to: ../Logs/3_filter_features_prevalence.log")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
