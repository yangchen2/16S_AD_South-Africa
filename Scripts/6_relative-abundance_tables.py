import os
import pandas as pd
from biom import load_table, Table
import logging
from biom.util import biom_open

# Setup logging
logging.basicConfig(
    filename='../Logs/6_relative-abundance_tables.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def biom_to_dataframe(biom_path):
    """
    Reads a BIOM table and converts it to a Pandas DataFrame.

    Parameters:
    biom_path (str): Path to the BIOM file.

    Returns:
    pd.DataFrame: DataFrame with the BIOM table contents.
    """
    try:
        logging.info(f"Reading BIOM table from {biom_path}.")
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe())
        logging.info(f"Successfully converted BIOM table from {biom_path} to DataFrame.")
        return df
    except Exception as e:
        logging.error(f"Error reading BIOM table at {biom_path}: {e}")
        return None

def convert_to_relative_abundance(df):
    """
    Converts a DataFrame to a relative abundance table.

    Parameters:
    df (pd.DataFrame): DataFrame to be converted.

    Returns:
    pd.DataFrame: DataFrame where each column is scaled to relative abundance.
    """
    try:
        logging.info("Converting DataFrame to relative abundance table.")
        # Remove columns with zero sum
        num_initial_columns = df.shape[1]
        df_nonzero = df.loc[:, df.sum(axis=0) != 0]
        num_removed_columns = num_initial_columns - df_nonzero.shape[1]
        logging.info(f"Removed {num_removed_columns} columns with zero sum.")

        # Convert to relative abundance
        df_relative_abundance = df_nonzero.div(df_nonzero.sum(axis=0), axis=1)
        logging.info("Successfully converted to relative abundance table.")
        logging.info(df_relative_abundance.head())
        return df_relative_abundance
    except Exception as e:
        logging.error(f"Error during conversion to relative abundance: {e}")
        return None

def save_dataframe_as_biom(df, output_path):
    """
    Saves a Pandas DataFrame as a BIOM table.

    Parameters:
    df (pd.DataFrame): DataFrame to save.
    output_path (str): Path to save the BIOM file.
    """
    try:
        logging.info(f"Saving DataFrame to BIOM format at {output_path}.")
        print(df.head())
        biom_table = Table(df.values, observation_ids=df.index, sample_ids=df.columns)
        with biom_open(output_path, 'w') as f:
            biom_table.to_hdf5(f, generated_by="Relative Abundance Script")
        logging.info(f"Successfully saved BIOM table to {output_path}.")
    except Exception as e:
        logging.error(f"Error saving DataFrame to BIOM format at {output_path}: {e}")

if __name__ == "__main__":
    # List of 7 BIOM table paths
    biom_table_paths = [
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Kingdom_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Phylum_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Class_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Order_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Family_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Genus_absolute.biom",
        "../Data/Tables/Absolute_Abundance_Tables/df_16S_filtered_feature_table_rare_Species_absolute.biom",
    ]

    output_dir = "../Data/Tables/Relative_Abundance_Tables"
    os.makedirs(output_dir, exist_ok=True)

    logging.info("Starting relative abundance table creation process.")

    for biom_path in biom_table_paths:
        try:
            logging.info(f"Processing {biom_path}...")
            # Step 1: Convert BIOM table to DataFrame
            df = biom_to_dataframe(biom_path)
            if df is not None:
                # Step 2: Convert to relative abundance table
                df_relative_abundance = convert_to_relative_abundance(df)

                # Save relative abundance table as BIOM
                output_file = os.path.join(output_dir, os.path.basename(biom_path).replace('_absolute.biom', '_relative_abundance.biom'))
                save_dataframe_as_biom(df_relative_abundance, output_file)
        except Exception as e:
            logging.error(f"Error processing {biom_path}: {e}")

    logging.info("Completed relative abundance table creation process.")
