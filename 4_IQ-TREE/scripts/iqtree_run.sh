#!/bin/bash

# ========== Check for argument ========== #
if [ -z "$1" ]; then
  echo "Usage: $0 <input_subfolder_name>"
  echo "Example: $0 filtered_50%Taxa_50Col_50%ColumnGapTolerance_50%RowGapTolerance"
  exit 1
fi

# ========== Configuration ========== #
# Project directory
PROJECT="/mnt/c/Users/beto2/Documents/Romalea_phylogenomics/4_IQ-TREE"

# Input folder (from command-line argument)
INPUT_SUBFOLDER="$1"

# Set directories
FILTER_DIR="$PROJECT/input/$INPUT_SUBFOLDER"
OUTPUT_DIR="$PROJECT/output/$INPUT_SUBFOLDER"
SCRIPT_DIR="$PROJECT/scripts"

# File names
ALIGNMENT_FILE="concatenated.fasta-out.nex"
PARTITION_FILE="partitions.txt"
MODIFIED_PARTITION_FILE="$FILTER_DIR/partitions_dna.txt"

# ========== Environment ========== #
source ~/miniconda3/etc/profile.d/conda.sh
conda activate iqtree_env  # Replace with your environment name

# ========== Setup ========== #
mkdir -p "$OUTPUT_DIR"

# Preprocess partition file: prepend "DNA, " to each line
sed 's/^/DNA, /' "$FILTER_DIR/$PARTITION_FILE" > "$MODIFIED_PARTITION_FILE"

# ========== IQ-TREE Analysis ========== #

# First: ModelFinder + Partition merging (no tree)
iqtree2 \
  -s "$FILTER_DIR/$ALIGNMENT_FILE" \
  -p "$MODIFIED_PARTITION_FILE" \
  -m MF+MERGE \
  -T AUTO \
  -pre "$OUTPUT_DIR/model_selection"

# Second: Ultrafast bootstrap with best scheme
iqtree2 \
  -s "$FILTER_DIR/$ALIGNMENT_FILE" \
  -p "$OUTPUT_DIR/model_selection.best_scheme.nex" \
  -B 1000 \
  -T AUTO \
  -pre "$OUTPUT_DIR/bootstrap_tree"
