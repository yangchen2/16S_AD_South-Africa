import pandas as pd
import qiime2 as q2
import biom
from biom.util import biom_open
from biom import load_table
import os
import glob
from collections import defaultdict
import logging

##########################################################################################
# SCRIPT 9: ASSIGN GENUS-ASV NAME TO EACH ASV FEATURE EXTRACTED FROM RF MODELS
##########################################################################################

# Setup logging
logging.basicConfig(
    filename='../logs/9_feature-name_to_rf_importance.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


# LOAD TAXONOMY DATA
logging.info("LOADING TAXONOMY DATA")
taxonomy_artifact = q2.Artifact.load('../Reference/2022.10.taxonomy.asv.tsv.qza')
gg_taxonomy = taxonomy_artifact.view(pd.DataFrame)

logging.info(f"Loaded taxonomy for {len(gg_taxonomy)} ASVs")
logging.info(f"Taxonomy columns: {list(gg_taxonomy.columns)}\n")


def create_asv_mapping_from_biom(tbl_path: str, level: str = 'Genus'):
    logging.info(f"Loading BIOM table from: {tbl_path}")
    biom_table = load_table(tbl_path)
    df = pd.DataFrame(biom_table.to_dataframe()).transpose()
    logging.info(f'Table shape: {df.shape}')

    df_tax = df.transpose()
    df_tax = df_tax.merge(gg_taxonomy['Taxon'], how='left', left_index=True, right_index=True)

    if df_tax['Taxon'].isnull().any():
        logging.warning(f"Some ASVs are missing taxonomy assignments at {level} level.")
        df_tax['Taxon'] = df_tax['Taxon'].fillna('Unknown')

    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_tax[taxonomy_levels] = df_tax['Taxon'].str.split(';', expand=True)
    df_tax = df_tax[df_tax[level].notnull()]
    logging.info(f"ASVs with {level} level taxonomy: {len(df_tax)}")

    sample_cols = [col for col in df_tax.columns if col not in taxonomy_levels + ['Taxon']]
    df_tax['TotalReadCount'] = df_tax[sample_cols].sum(axis=1)

    index_map = {}
    abundance_data = []

    for taxon, group in df_tax.groupby(level):
        group_sorted = group.sort_values('TotalReadCount', ascending=False)
        for i, idx in enumerate(group_sorted.index, 1):
            clean_taxon = taxon.strip().replace(' ', '_')
            new_id = f"{clean_taxon}_ASV-{i}"
            index_map[idx] = new_id
            abundance_data.append({
                'ASV_Sequence': idx,
                'ASV_Name': new_id,
                'Genus': taxon,
                'TotalReadCount': group_sorted.loc[idx, 'TotalReadCount'],
                'Full_Taxonomy': df_tax.loc[idx, 'Taxon']
            })

    df_with_abundance = pd.DataFrame(abundance_data)
    df_with_abundance = df_with_abundance.sort_values('TotalReadCount', ascending=False)

    logging.info(f"Created mapping for {len(index_map)} ASVs")
    logging.info(f"Number of unique {level}: {df_tax[level].nunique()}")

    return index_map, df_with_abundance


logging.info("CREATING ASV MAPPING FROM BIOM TABLE")

biom_table_path = "../Data/Tables/Count_Tables/2_209766_feature_table_dedup.biom"

if 'df' in globals() and not os.path.exists(biom_table_path):
    logging.info("\nUsing DataFrame 'df' to create mapping...")
    asv_abundances = df.sum(axis=0).sort_values(ascending=False)
    index_map = {}
    abundance_data = []
    df_tax = gg_taxonomy.copy()
    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_tax[taxonomy_levels] = df_tax['Taxon'].str.split(';', expand=True)
    for genus, group in df_tax.groupby('Genus'):
        genus_asvs = [asv for asv in group.index if asv in asv_abundances.index]
        if not genus_asvs:
            continue
        genus_abundances = asv_abundances[genus_asvs].sort_values(ascending=False)
        for i, (asv_seq, abundance) in enumerate(genus_abundances.items(), 1):
            clean_genus = genus.strip().replace(' ', '_')
            new_id = f"{clean_genus}_ASV-{i}"
            index_map[asv_seq] = new_id
            abundance_data.append({
                'ASV_Sequence': asv_seq,
                'ASV_Name': new_id,
                'Genus': genus,
                'TotalReadCount': abundance,
                'Full_Taxonomy': df_tax.loc[asv_seq, 'Taxon']
            })
    df_with_abundance = pd.DataFrame(abundance_data)
    df_with_abundance = df_with_abundance.sort_values('TotalReadCount', ascending=False)
    logging.info(f"Created mapping for {len(index_map)} ASVs from DataFrame")

else:
    index_map, df_with_abundance = create_asv_mapping_from_biom(biom_table_path, level='Genus')

logging.info(f"\nTotal ASVs mapped: {len(index_map)}")
logging.info(f"Total unique genera: {df_with_abundance['Genus'].nunique()}")
logging.info(df_with_abundance[['ASV_Name', 'Genus', 'TotalReadCount']].head(25).to_string(index=False))


def add_readable_names_to_csv(csv_path, asv_mapping, output_path=None):
    df = pd.read_csv(csv_path, index_col=0)

    # Sort by mean_importance if present
    if 'mean_importance' in df.columns:
        df = df.sort_values('mean_importance', ascending=False)

    # Add readable names
    df_updated = df.copy()
    df_updated.insert(0, 'ASV_Name', df_updated.index.map(asv_mapping))

    # Remove rows where ASV_Name starts with 'g___ASV' or is NA
    df_updated = df_updated[
        df_updated['ASV_Name'].notna() &
        ~df_updated['ASV_Name'].str.startswith('g___ASV', na=False)
    ]

    mapped_count = df_updated['ASV_Name'].notna().sum()
    total_count = len(df_updated)

    if output_path is None:
        base_name = os.path.splitext(csv_path)[0]
        output_path = f"{base_name}_ASV-name_known.csv"

    # Overwrite if file exists
    df_updated.to_csv(output_path, index=True)

    return df_updated, mapped_count, total_count


logging.info("PROCESSING FEATURE IMPORTANCE CSV FILES")

fi_directory = '../Data/RF_Feature_Importances'
os.makedirs(fi_directory, exist_ok=True)

csv_files = glob.glob(os.path.join(fi_directory, '*.csv'))
csv_files = [f for f in csv_files if 'readable' not in f.lower() 
             and '_named' not in f.lower() 
             and 'mapping_summary' not in f.lower()]

logging.info(f"\nFound {len(csv_files)} CSV files to process:")
for csv_file in csv_files:
    logging.info(f"  - {os.path.basename(csv_file)}")

if len(csv_files) == 0:
    logging.info("\nNo CSV files found to process.")
else:
    logging.info("PROCESSING FILES")
    summary_results = []

    for csv_file in csv_files:
        filename = os.path.basename(csv_file)
        logging.info(f"\n{filename}")

        try:
            df_updated, mapped_count, total_count = add_readable_names_to_csv(
                csv_file, index_map
            )
            logging.info(f"  Mapped {mapped_count}/{total_count} ASVs (excluding g__ASV entries)")

            base_name = os.path.splitext(filename)[0]
            output_filename = f"{base_name}_ASV-name.csv"
            logging.info(f"  Saved to: {output_filename}")

            summary_results.append({
                'Original_File': filename,
                'Output_File': output_filename,
                'Total_ASVs': total_count,
                'Mapped_ASVs': mapped_count,
                'Unmapped_ASVs': total_count - mapped_count,
                'Mapping_Rate': f"{mapped_count/total_count*100:.1f}%"
            })

            logging.info("  Top 5 most important ASVs:")
            for idx, row in df_updated.head(5).iterrows():
                readable_name = row['ASV_Name']
                importance = row.get('mean_importance', 'N/A')
                if idx in df_with_abundance['ASV_Sequence'].values:
                    abundance = df_with_abundance[df_with_abundance['ASV_Sequence']==idx]['TotalReadCount'].values[0]
                    logging.info(f"    {readable_name}: importance={importance:.6f}, abundance={abundance:.6f}")
                else:
                    logging.info(f"    {readable_name}: importance={importance}")

        except Exception as e:
            logging.info(f"  ERROR: {str(e)}")
            logging.error(f"Failed to process {filename}: {str(e)}")
            summary_results.append({
                'Original_File': filename,
                'Output_File': 'FAILED',
                'Total_ASVs': 0,
                'Mapped_ASVs': 0,
                'Unmapped_ASVs': 0,
                'Mapping_Rate': '0%'
            })

os.makedirs('../Data/Taxonomy', exist_ok=True)
mapping_path = '../Data/Taxonomy/ASV_readable_name_mapping_abundance_ranked.csv'
df_with_abundance.to_csv(mapping_path, index=False)
logging.info(f"Saved detailed mapping: {mapping_path}")
logging.info("COMPLETE!")
