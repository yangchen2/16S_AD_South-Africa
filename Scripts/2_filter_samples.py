import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import logging
import os

##########################################################################################
# SCRIPT 2: FILTERS RAW BIOM TABLE TO MATCH DEDUPLICATED SAMPLE SET IN SUBSET METADATA
##########################################################################################

# Setup logging
os.makedirs('../Logs', exist_ok=True)
logging.basicConfig(
    filename='../Logs/2_filter-samples.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def read_and_convert_biom(biom_path: str) -> pd.DataFrame:
    """
    Load BIOM file and convert to pandas DataFrame (samples as rows, features as columns)
    """
    try:
        logging.info(f"Loading BIOM file: {biom_path}")
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe()).T
        df.index = df.index.astype(str)
        logging.info(f"BIOM table shape: {df.shape}")
        return df
    except Exception as e:
        logging.error(f"Error in loading BIOM: {e}")
        raise

def load_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Load sample metadata and set #sample-id as index, harmonizing underscores
    """
    try:
        metadata = pd.read_csv(metadata_path, sep='\t')
        metadata['#sample-id'] = metadata['#sample-id'].astype(str).str.replace('_', '')
        metadata.set_index('#sample-id', inplace=True)
        logging.info(f"Metadata loaded with shape: {metadata.shape}")
        return metadata
    except Exception as e:
        logging.error(f"Error in loading metadata: {e}")
        raise

def subset_biom_metadata(df: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Subset BIOM table to only include samples present in metadata
    """
    logging.info(f"Original BIOM table samples: {len(df)}")
    logging.info(f"Metadata samples: {len(metadata)}")

    # Remove numeric prefix (e.g., "15564.") if present
    df.index = df.index.str.replace(r'^15564\.', '', regex=True)

    # Subset to overlap
    df_subset = df[df.index.isin(metadata.index)]
    logging.info(f"Shape after metadata subset: {df_subset.shape}")
    logging.info(f"Samples removed: {len(df) - len(df_subset)}")

    return df_subset

def save_as_biom(df: pd.DataFrame, output_path: str):
    """
    Save DataFrame as BIOM format file
    """
    df = df.transpose()
    table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with biom_open(output_path, 'w') as f:
        table.to_hdf5(f, "filtered sample script")
    logging.info(f"Saved filtered BIOM table to {output_path}")

def main():
    logging.info("Starting sample filtering script")

    biom_path = '../Data/Tables/Count_Tables/1_209766_feature_table.biom'
    metadata_path = '../Metadata/16S_AD_South-Africa_metadata_subset.tsv'
    output_path = '../Data/Tables/Count_Tables/2_209766_feature_table_dedup.biom'

    try:
        # Load BIOM table and metadata
        df = read_and_convert_biom(biom_path)
        metadata = load_metadata(metadata_path)

        # Subset BIOM table to match metadata samples
        df_filtered = subset_biom_metadata(df, metadata)

        # Save filtered BIOM
        save_as_biom(df_filtered, output_path)

        logging.info("Sample filtering completed successfully")

    except Exception as e:
        logging.error(f"Script failed with error: {e}")
        raise

if __name__ == "__main__":
    main()
