import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
from numpy.random import RandomState
import matplotlib.pyplot as plt
import logging

# Setup logging
logging.basicConfig(filename='../Logs/2_filter-samples_and_rarefy.log', level=logging.INFO,
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
        df = pd.DataFrame(biom_table.to_dataframe()).T
        df.index = df.index.str.replace('^15564\.', '', regex=True)  # Remove '15564.' prefix from column names
        logging.info(f"BIOM table shape: {df.shape}")
        logging.info(df.head())
        logging.info("BIOM file successfully converted to DataFrame")
        
        return df
    except Exception as e:
        logging.error(f"Error in processing BIOM file: {e}")
        raise

def load_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Load metadata and set up for merging with the BIOM data.

    Parameters:
    metadata_path (str): Path to the metadata file.

    Returns:
    pd.DataFrame: Metadata with '#sample-id' as index.
    """
    try:
        metadata = pd.read_csv(metadata_path, sep='\t')
        metadata.set_index('#sample-id', inplace=True)
        logging.info(f"Metadata loaded with shape: {metadata.shape}")
        return metadata
    except Exception as e:
        logging.error(f"Error in loading metadata: {e}")
        raise



def plot_sorted_read_counts(biom_df: pd.DataFrame, metadata: pd.DataFrame):
    """
    Plot the read counts per sample in descending order, colored by specific case types.

    Parameters:
    biom_df (pd.DataFrame): Input DataFrame representing the BIOM table.
    metadata (pd.DataFrame): DataFrame containing metadata for each sample.
    """
    # Calculate the total read counts per sample
    sample_read_counts = biom_df.sum(axis=1).sort_values(ascending=False)
    
    # Merge read counts with metadata and drop NaN case_types
    sample_read_counts_df = sample_read_counts.to_frame(name='read_count').merge(metadata[['case_type']], 
                                                                                 left_index=True, right_index=True, 
                                                                                 how='left')
    sample_read_counts_df = sample_read_counts_df.dropna(subset=['case_type'])  # Remove rows where case_type is NaN

    # Assign a numeric index for plotting in order
    sample_read_counts_df = sample_read_counts_df.reset_index(drop=True)
    
    # Define specific colors for each case type
    color_map = {
        'control-nonlesional skin': 'blue',
        'case-nonlesional skin': 'orange',
        'case-lesional skin': 'red',
        'control-anterior nares': 'green',
        'case-anterior nares': 'pink'
    }
    sample_read_counts_df['color'] = sample_read_counts_df['case_type'].map(color_map)

    # Plot read counts with colors by case type
    plt.figure(figsize=(10, 6))
    plt.scatter(x=sample_read_counts_df.index, 
                y=sample_read_counts_df['read_count'], 
                c=sample_read_counts_df['color'], 
                s=20)

    # Add legend for case types in the specified order
    legend_order = ['control-nonlesional skin', 'case-nonlesional skin', 'case-lesional skin', 
                    'control-anterior nares', 'case-anterior nares']
    for case_type in legend_order:
        plt.scatter([], [], color=color_map[case_type], label=case_type)
    
    plt.title('Read Counts per Sample (Sorted High to Low)')
    plt.xlabel('Sample (sorted by read count)')
    plt.ylabel('Read Count')
    plt.xticks([])  # Remove x-axis tick labels
    plt.legend(title="Case Type")
    plt.grid(True)
    plt.savefig('../Plots/read_counts_by_case_type.png', dpi=600)
    logging.info("Read counts per sample plotted by case type in descending order.")





if __name__ == '__main__':
    try:
        # Paths for rs210 per-genome table and metadata
        biom_path = '../Data/Tables/209766_filtered_feature_table.biom'
        metadata_path = '../Data/Metadata/updated_clean_ant_skin_metadata.tab'
        output_path = '../Data/Tables/209766_filtered_feature_table_rare.biom'

        # Read and convert BIOM to DataFrame
        df = read_and_convert_biom(biom_path)
        
        # Load metadata
        metadata = load_metadata(metadata_path)
        
        # Plot read counts per sample to help decide on filtering threshold, colored by sample type
        plot_sorted_read_counts(df, metadata)
        
        # Further processing steps such as rarefaction can go here as per original code

        logging.info("Script completed successfully.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
