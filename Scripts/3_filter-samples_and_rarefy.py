import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
from numpy.random import RandomState
import matplotlib.pyplot as plt
import logging

# Setup logging
logging.basicConfig(filename='../Logs/3_filter-samples_and_rarefy.log', level=logging.INFO,
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
    plt.savefig('../Plots/Dataset_stats/read_counts_by_case_type.png', dpi=600)
    logging.info("Read counts per sample plotted by case type in descending order.")



def plot_rarefaction_depths(biom_df: pd.DataFrame) -> int:
    """
    Plot the distribution of sample read depths to visualize possible rarefying depths
    and return the mean depth.

    Parameters:
    biom_df (pd.DataFrame): Input DataFrame representing the BIOM table.

    Returns:
    int: Mean depth across samples.
    """
    # Calculate sample read depths
    sample_depths = biom_df.sum(axis=1)
    
    # Calculate mean depth
    mean_depth = int(sample_depths.mean())
    
    # Plot histogram of read depths
    plt.figure(figsize=(10, 6))
    plt.hist(sample_depths, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(sample_depths.min(), color='red', linestyle='--', label=f'Min Depth: {sample_depths.min()}')
    plt.axvline(mean_depth, color='blue', linestyle='--', label=f'Mean Depth: {mean_depth}')
    plt.title('Distribution of Sample Read Depths')
    plt.xlabel('Read Depth')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig('../Plots/Dataset_stats/read_counts_by_sample.png', dpi=600)
    logging.info(f"Rarefaction depth distribution plotted. Min depth: {sample_depths.min()}. Mean depth: {mean_depth}")
    return mean_depth



def rarefy_table(biom_df: pd.DataFrame, seed: int = 42, depth: int = None) -> pd.DataFrame:
    """
    Rarefy a single BIOM table to a specified sampling depth.

    Parameters:
    biom_df (pd.DataFrame): Input DataFrame representing the BIOM table.
    seed (int): Seed for the random number generator to ensure reproducibility.
    depth (int): Depth to which to rarefy samples. Defaults to minimum sample depth if not provided.

    Returns:
    pd.DataFrame: A rarefied DataFrame.
    """
    # Calculate minimum sampling depth if no custom depth is provided
    if depth is None:
        depth = biom_df.sum(axis=1).min().astype(int)
    logging.info(f"Using rarefaction depth: {depth}")

    # Rarefaction process
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



# def filter_low_depth_samples(biom_df: pd.DataFrame, min_depth: int = 1000) -> pd.DataFrame:
#     """
#     Filter out samples with a read depth below a specified threshold.

#     Parameters:
#     biom_df (pd.DataFrame): Input DataFrame representing the BIOM table.
#     min_depth (int): Minimum read depth threshold.

#     Returns:
#     pd.DataFrame: Filtered DataFrame with only samples having read depth above the threshold.
#     """
#     initial_sample_count = biom_df.shape[0]
#     filtered_df = biom_df[biom_df.sum(axis=1) >= min_depth]
#     final_sample_count = filtered_df.shape[0]

#     # Log the number of samples before and after filtering
#     logging.info(f"Filtering samples with read depth below {min_depth}")
#     logging.info(f"Number of samples before filtering: {initial_sample_count}")
#     logging.info(f"Number of samples after filtering: {final_sample_count}")
#     logging.info(f"Number of samples removed: {initial_sample_count - final_sample_count}")

#     return filtered_df



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
        table.to_hdf5(f, "rarefaction script")
    logging.info(f"DataFrame saved as BIOM table to {output_path}")



if __name__ == '__main__':
    try:
        # Paths for rs210 per-genome table
        biom_path = '../Data/Tables/209766_filtered_feature_table.biom'
        metadata_path = '../Data/Metadata/updated_clean_ant_skin_metadata.tab'
        output_path = '../Data/Tables/209766_filtered_feature_table_rare.biom'

        # Read and convert BIOM to DataFrame
        df = read_and_convert_biom(biom_path)
        
        # Load metadata
        metadata = load_metadata(metadata_path)

        # Filter out samples with low read depth
        # df = filter_low_depth_samples(df, min_depth=500)

        # Plot read counts per sample to help decide on filtering threshold
        plot_sorted_read_counts(df, metadata)
        
        # Plot rarefaction depths to choose rarefying depth
        mean_depth = plot_rarefaction_depths(df)

        # Perform rarefaction using mean depth
        rarefied_df = rarefy_table(df, seed=42, depth=350)  # based on qiime2 rarefaction curve
        
        # Save the rarefied table as a BIOM file
        save_as_biom(rarefied_df, output_path)
        
        logging.info("Script completed successfully.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

