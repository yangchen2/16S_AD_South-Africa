import os
import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import scipy.stats as ss
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2 as q2
import warnings
import logging

warnings.simplefilter(action='ignore', category=FutureWarning)

# Setup logging
logging.basicConfig(
    filename='../logs/5_taxonomy_tbl_collapse.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Read in greengenes2 taxonomy mapping file
logging.info("Loading GreenGreens2 taxonomy... (this will take ~2 minutes)")
gg_taxonomy = q2.Artifact.load('../Reference/2022.10.taxonomy.asv.tsv.qza').view(pd.DataFrame)
logging.info("GreenGreens2 taxonomy loaded")


def group_by_all_taxonomy_levels(tbl_path: str):
    """
    Read rarefied 16S, attach gg2 taxonomy info, and collapse by taxonomy level to multiple BIOMs.

    Parameters:
    tbl_path (str): Path to the BIOM table

    Returns:
    dict: A dictionary containing 7 dataframes for 7 taxonomy levels
    """
    logging.info("STEP 1: Converting BIOMs to dataframes.")
    # Read in BIOM tables and convert to Pandas df
    biom_table_16S = load_table(tbl_path)
    df_16S = pd.DataFrame(biom_table_16S.to_dataframe()).transpose()

    # Dictionary to hold the grouped dataframes for each table and each level
    grouped_dfs = {'df_16S': {}}

    logging.info("STEP 2: Attaching taxonomy info and collapsing by taxonomy level.")
    os.makedirs("../Intermediate_Outputs", exist_ok=True)

    # Process each DataFrame
    for df_name, df in zip(['df_16S'], [df_16S]):
        # Transpose df
        df = df.transpose()

        # Fill in NA values
        df.fillna(0, inplace=True)

        # Cast to int
        df = df.astype(int)

        # Merge taxonomy column based on common index
        unmatched_features = set(df.index) - set(gg_taxonomy.index)
        if unmatched_features:
            logging.warning(f"{len(unmatched_features)} features in the BIOM table do not have corresponding taxonomy information.")
        df = df.merge(gg_taxonomy['Taxon'], how='left', left_index=True, right_index=True)
        df.to_csv(f"../Intermediate_Outputs/{df_name}_with_taxonomy.csv")

        # Define all taxonomy levels
        taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

        logging.info("Looping through each taxonomy level and performing grouping.")
        for level in taxonomy_levels:
            try:
                # Copy the original df to avoid modifying it permanently
                df_level = df.copy()

                # Split Taxon column
                if 'Taxon' not in df_level.columns or df_level['Taxon'].isnull().any():
                    logging.warning(f"Some rows are missing taxonomy information for level {level}. These will be excluded.")
                
                df_level[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = df_level['Taxon'].str.split(';', expand=True)

                # Drop unnecessary taxonomy levels
                levels_to_remove = taxonomy_levels.copy()
                levels_to_remove.remove(level)
                df_reduced = df_level.drop(columns=levels_to_remove)

                # Group by (collapse) features to the current taxonomy level
                df_grouped = df_reduced.groupby(level).sum(numeric_only=True)

                # Sort by row sums
                df_grouped = df_grouped.loc[df_grouped.sum(axis=1).sort_values(ascending=False).index]

                # Save BIOM
                obs_ids = df_grouped.index
                samp_ids = df_grouped.columns
                biom_table = biom.table.Table(df_grouped.values, observation_ids=obs_ids, sample_ids=samp_ids)
                print('Table shape: ' + str(biom_table.shape))
                biom_output_file = f"../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_1pct_rare_{level}.biom"
                os.makedirs("../Data/Tables/Absolute_Abundance_Tables", exist_ok=True)
                with biom_open(biom_output_file, 'w') as f:
                    biom_table.to_hdf5(f, generated_by="Collapsed by taxonomy level")

                logging.info(f"Saved collapsed BIOM at {level} level: {biom_output_file}")

            except Exception as e:
                logging.error(f"Error processing taxonomy level {level}: {e}")


    logging.info("-----> COMPLETED: Taxonomy collapsed BIOM files successfully converted to DataFrames")
    return grouped_dfs


if __name__ == '__main__':
    try:
        # File paths for the BIOM files
        biom_16S_path = "../Data/Tables/Absolute_Abundance_Tables/209766_filtered_by_prevalence_1pct_rare.biom"

        # Create collapsed taxa BIOMs
        grouped_dfs = group_by_all_taxonomy_levels(biom_16S_path)
    except Exception as e:
        logging.error(f"An error occurred in the main execution: {e}")
