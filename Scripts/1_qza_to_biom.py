import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import logging

# Set up logging
logging.basicConfig(filename='../Logs/1_qza_to_biom.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def qza_to_biom(qza_path: str):
    try:
        # Convert from qza file to pandas df
        biom_df = q2.Artifact.load(qza_path).view(pd.DataFrame).transpose()

        # Check df contents
        logging.info(f"Table shape: {biom_df.shape}")
        logging.info(biom_df.head())

        # Save as BIOM file
        obs_ids = biom_df.index
        samp_ids = biom_df.columns

        biom_table = biom.table.Table(biom_df.values, observation_ids=obs_ids, sample_ids=samp_ids)
        
        # Remove the .qza extension and append .biom
        biom_output_file = qza_path.rsplit('.qza', 1)[0] + '.biom'

        with biom_open(biom_output_file, 'w') as f:
            biom_table.to_hdf5(f, generated_by="qza_to_biom.py")
        
        # Log the success
        logging.info(f"{qza_path} successfully converted to {biom_output_file}")
        
    except Exception as e:
        # Log any errors
        logging.error(f"An error occurred while converting {qza_path}: {e}")
        raise



if __name__ == '__main__':
    try:
        # Path for rs210 per-genome table
        qza_path = '../Data/Tables/209766_filtered_feature_table.qza'

        # Convert
        qza_to_biom(qza_path)

    except Exception as e:
        logging.error(f"An error occurred: {e}")