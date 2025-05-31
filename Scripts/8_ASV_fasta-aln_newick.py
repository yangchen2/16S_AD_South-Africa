#!/usr/bin/env python

import subprocess
import sys
import logging
import os

# Setup logging
os.makedirs("../Logs", exist_ok=True)
logging.basicConfig(
    filename='../Logs/8_ASV_fasta-aln_newick.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define prevalence thresholds to process
prevalence_labels = ["1pct", "5pct", "10pct"]

# Define file path templates
fasta_dir = "../Data/Fasta"
tree_dir = "../Data/Trees"

def run_mafft(input_fasta, output_aln):
    """Run MAFFT to generate multiple sequence alignment."""
    print(f"Running MAFFT on {input_fasta}...")
    logging.info(f"Running MAFFT on {input_fasta}...")
    try:
        with open(output_aln, 'w') as fout:
            subprocess.run(["mafft", "--auto", input_fasta], stdout=fout, check=True)
        print(f"Alignment written to {output_aln}")
        logging.info(f"Alignment written to {output_aln}")
    except subprocess.CalledProcessError:
        print("Error: MAFFT alignment failed.")
        logging.error("MAFFT alignment failed.")
        sys.exit(1)

def run_fasttree(input_aln, output_tree):
    """Run FastTree to build a phylogenetic tree."""
    print(f"Running FastTree on {input_aln}...")
    logging.info(f"Running FastTree on {input_aln}...")
    try:
        with open(output_tree, 'w') as fout:
            subprocess.run(["FastTree", "-nt", input_aln], stdout=fout, check=True)
        print(f"Tree written to {output_tree}")
        logging.info(f"Tree written to {output_tree}")
    except subprocess.CalledProcessError:
        print("Error: FastTree failed.")
        logging.error("FastTree failed.")
        sys.exit(1)

if __name__ == '__main__':
    for label in prevalence_labels:
        input_fasta = f"{fasta_dir}/209766_filtered_by_prevalence_{label}_Healthy_only.fasta"
        output_aln = f"{fasta_dir}/209766_filtered_by_prevalence_{label}_Healthy_only.aln.fasta"
        output_tree = f"{tree_dir}/209766_filtered_by_prevalence_{label}_Healthy_only.tree.nwk"

        if os.path.exists(input_fasta):
            run_mafft(input_fasta, output_aln)
            run_fasttree(output_aln, output_tree)
        else:
            print(f"Input file {input_fasta} not found. Skipping.")
            logging.warning(f"Input file {input_fasta} not found. Skipping.")
