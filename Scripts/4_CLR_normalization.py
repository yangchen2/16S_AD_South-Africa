import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
import logging

# Setup logging
logging.basicConfig(filename='../Logs/4_clr_normalization.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def read_and_convert_biom(biom_path: str) -> pd.DataFrame:
    """
    Read a BIOM-format table and convert it into a Pandas dataframe.

    Parameters:
    biom_path (str): The file path to the biom table.

    Returns:
    pd.DataFrame: The dataframe corresponding to the biom table.
    """
    logging.info(f"Loading BIOM file: {biom_path}")
    try:
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe().T)
        logging.info(f"BIOM table shape: {df.shape}")
        logging.info("BIOM file successfully converted to DataFrame")
        return df
    except Exception as e:
        logging.error(f"Error in processing BIOM file: {e}")
        raise



def clr_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform Centered Log-Ratio (CLR) transformation on a DataFrame.
    
    Parameters:
    df (pd.DataFrame): DataFrame where rows represent samples and columns represent features.
    
    Returns:
    pd.DataFrame: CLR-transformed DataFrame.
    """
    # Add a pseudocount to avoid zeros
    df += 1
    
    # Compute geometric mean for each sample (row-wise)
    geometric_means = df.apply(lambda x: np.exp(np.log(x).mean()), axis=1)
    
    # Perform CLR transformation
    clr_transformed_df = np.log(df.div(geometric_means, axis=0))
    logging.info(clr_transformed_df.head())
    
    logging.info("CLR transformation completed.")
    return clr_transformed_df



def save_as_biom(df: pd.DataFrame, output_path: str):
    """
    Save a pandas DataFrame as a BIOM table.

    Parameters:
    df (pd.DataFrame): The DataFrame to save.
    output_path (str): Path to the output BIOM file.
    """
    df = df.transpose()
    table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    with biom_open(output_path, 'w') as f:
        table.to_hdf5(f, "clr transformation script")
    logging.info(f"DataFrame saved as BIOM table to {output_path}")

if __name__ == '__main__':
    try:
        # Paths for input and output files
        biom_path = '../Data/Tables/209766_filtered_feature_table.biom'
        output_path = '../Data/Tables/209766_filtered_feature_table_clr.biom'

        # Read and convert BIOM to DataFrame
        df = read_and_convert_biom(biom_path)
        
        # Perform CLR normalization
        clr_df = clr_transform(df)
        
        # Save the CLR-normalized table as a BIOM file
        save_as_biom(clr_df, output_path)
        
        logging.info("Script completed successfully.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
