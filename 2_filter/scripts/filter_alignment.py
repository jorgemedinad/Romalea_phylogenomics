#!/usr/bin/env python3
"""
Alignment Filter Script (with CLI arguments)
Author: Jorge H. Medina-Duran
Converted to Python: May 2025

This script filters DNA alignments in FASTA format based on:
- Maximum allowed gaps per column
- Maximum allowed gaps per taxon (row)
- Minimum number of alignment columns
- Minimum number of taxa retained

All parameters are passed via command-line arguments.

Directory Structure:
  working_directory/
  ├── input/   # Contains unfiltered .fasta alignments
  ├── output/  # Filtered alignments will be saved here
  └── scripts/ # This script should be located here
"""

import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np

# ----- Parse CLI arguments -----
parser = argparse.ArgumentParser(description="Filter alignments based on gap thresholds and taxon coverage.")
parser.add_argument("--total_taxa", type=int, required=True, help="Total number of taxa in the alignment.")
parser.add_argument("--column_gap", type=float, default=0.25, help="Max proportion of gaps allowed per column.")
parser.add_argument("--row_gap", type=float, default=0.25, help="Max proportion of gaps allowed per row (taxon).")
parser.add_argument("--min_columns", type=int, default=50, help="Minimum number of columns after filtering.")
parser.add_argument("--taxa_percent", type=float, default=0.5, help="Minimum proportion of taxa to retain.")

args = parser.parse_args()

# ----- Filtering parameters -----
total_taxa = args.total_taxa
column_gap_threshold = args.column_gap
taxa_gap_threshold = args.row_gap
min_columns = args.min_columns
taxa_percent = args.taxa_percent
min_taxa = int(total_taxa * taxa_percent)

# ----- Set working directories -----
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, ".."))
input_dir = os.path.join(base_dir, "input")
output_base = os.path.join(base_dir, "output")

output_folder_name = f"filtered_{int(taxa_percent*100)}%Taxa_{min_columns}Col_"                      f"{int(column_gap_threshold*100)}%ColumnGapTolerance_"                      f"{int(taxa_gap_threshold*100)}%RowGapTolerance"
final_output_dir = os.path.join(output_base, output_folder_name)

os.makedirs(final_output_dir, exist_ok=True)

# ----- Process each FASTA file -----
for filename in os.listdir(input_dir):
    if not filename.endswith(".fasta"):
        continue

    input_path = os.path.join(input_dir, filename)
    records = list(SeqIO.parse(input_path, "fasta"))

    if len(records) == 0:
        print(f"Empty or invalid FASTA file: {filename}")
        continue

    taxa_names = [rec.id for rec in records]
    seq_array = np.array([list(str(rec.seq)) for rec in records])
    num_taxa, num_columns = seq_array.shape

    # Step 1: Filter columns by gap threshold
    gap_counts = np.sum((seq_array == "-") | (seq_array == "?"), axis=0)
    columns_to_keep = gap_counts < (column_gap_threshold * num_taxa)
    
    if np.sum(columns_to_keep) < min_columns:
        print(f"File {filename} has fewer than {min_columns} columns after filtering. Skipping.")
        continue

    seq_array = seq_array[:, columns_to_keep]
    num_columns = seq_array.shape[1]

    # Step 2: Filter rows (taxa) by gap threshold
    gap_counts_per_row = np.sum((seq_array == "-") | (seq_array == "?"), axis=1)
    taxa_to_keep = gap_counts_per_row < (taxa_gap_threshold * num_columns)

    if np.sum(taxa_to_keep) == 0:
        print(f"No taxa meet the gap threshold in file {filename}. Skipping.")
        continue

    seq_array = seq_array[taxa_to_keep, :]
    final_taxa_names = np.array(taxa_names)[taxa_to_keep]

    # Step 3: Ensure minimum number of taxa
    if len(final_taxa_names) < min_taxa:
        print(f"File {filename} has fewer than {min_taxa} taxa after filtering. Skipping.")
        continue

    # Step 4: Write output
    filtered_records = []
    for i, name in enumerate(final_taxa_names):
        seq_str = "".join(seq_array[i])
        filtered_records.append(SeqRecord(seq=Seq(seq_str), id=name, description=""))

    output_path = os.path.join(final_output_dir, filename.replace(".fasta", "_filtered_final.fasta"))
    SeqIO.write(filtered_records, output_path, "fasta")
    print(f"Filtered alignment saved to {output_path}")
