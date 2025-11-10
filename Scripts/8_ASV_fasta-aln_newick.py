#!/usr/bin/env python

import subprocess
import sys
import logging
import os
from Bio import SeqIO

##########################################################################################
# SCRIPT 8: ALIGN ASVs IN FASTA AND CREATE NEWICK TREES (ALL, SKIN, NASAL)
##########################################################################################

# Setup logging
os.makedirs("../Logs", exist_ok=True)
logging.basicConfig(
    filename="../Logs/8_ASV_fasta-aln_newick.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Directories
fasta_dir = "../Data/Fasta"
tree_dir = "../Data/Trees"
os.makedirs(tree_dir, exist_ok=True)


def run_mafft(input_fasta, output_aln):
    """Run MAFFT to generate multiple sequence alignment."""
    logging.info(f"Running MAFFT on {input_fasta}")
    print(f"Running MAFFT on {input_fasta}...")
    try:
        with open(output_aln, "w") as fout:
            subprocess.run(["mafft", "--auto", input_fasta],
                           stdout=fout, stderr=subprocess.PIPE, check=True)
        print(f"Alignment written to {output_aln}")
        logging.info(f"Alignment successfully written to {output_aln}")
    except subprocess.CalledProcessError as e:
        print(f"MAFFT failed on {input_fasta}: {e}")
        logging.error(f"MAFFT alignment failed on {input_fasta}: {e}")


def run_fasttree(input_aln, output_tree):
    """Run FastTree to build a phylogenetic tree."""
    logging.info(f"Running FastTree on {input_aln}")
    print(f"Running FastTree on {input_aln}...")
    try:
        with open(output_tree, "w") as fout:
            subprocess.run(["FastTree", "-nt", input_aln],
                           stdout=fout, stderr=subprocess.PIPE, check=True)
        print(f"Tree written to {output_tree}")
        logging.info(f"Tree successfully written to {output_tree}")
    except subprocess.CalledProcessError as e:
        print(f"FastTree failed on {input_aln}: {e}")
        logging.error(f"FastTree failed on {input_aln}: {e}")


def count_sequences(fasta_path):
    """Count sequences in FASTA."""
    try:
        return sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
    except Exception:
        return 0


if __name__ == "__main__":
    prevalence_thresholds = ["10pct", "5pct", "1pct", '0pct']
    rarefaction_depths = [350, 1000, 1500, 2000]  # Based on QIIME2 rarefaction curve
    specimen_variants = ["all", "skin", "nasal"]

    logging.info("Starting MAFFT + FastTree for ASV FASTAs (ALL, SKIN, NASAL).")

    for threshold in prevalence_thresholds:
        for depth in rarefaction_depths:
            print(f"\n=== Processing prevalence {threshold}, depth {depth} ===")

            for variant in specimen_variants:
                input_fasta = (
                    f"{fasta_dir}/209766_feature_table_dedup_prev-filt-"
                    f"{threshold}_rare-{depth}_Genus-ASV_{variant}.fasta"
                )
                output_aln = (
                    f"{fasta_dir}/209766_feature_table_dedup_prev-filt-"
                    f"{threshold}_rare-{depth}_Genus-ASV_{variant}_aln.fasta"
                )
                output_tree = (
                    f"{tree_dir}/209766_feature_table_dedup_prev-filt-"
                    f"{threshold}_rare-{depth}_Genus-ASV_{variant}_aln.nwk"
                )

                # Skip missing FASTAs
                if not os.path.exists(input_fasta):
                    print(f"FASTA not found: {input_fasta}. Skipping.")
                    logging.warning(f"FASTA not found for {variant} ({threshold}, depth={depth}). Skipping.")
                    continue

                # Check sequence count
                seq_count = count_sequences(input_fasta)
                if seq_count < 2:
                    print(f"Skipping {variant} — only {seq_count} sequences found.")
                    logging.warning(f"Skipping {variant} ({threshold}, depth={depth}) — insufficient sequences ({seq_count}).")
                    continue

                # Run MAFFT + FastTree
                run_mafft(input_fasta, output_aln)
                run_fasttree(output_aln, output_tree)

                logging.info(f"Completed tree for {variant} ({threshold}, depth={depth})")

    logging.info("All alignments and trees generated successfully (ALL, SKIN, NASAL).")
    print("\nAll alignments and trees generated successfully.")
