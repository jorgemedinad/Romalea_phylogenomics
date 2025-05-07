#!/bin/bash

#project directory
PROJECT="/mnt/c/Users/beto2/Documents/Romalea_phylogenomics/pipeline/3_AMAS"  ### Change this to your project directory

# Input directory (full path to filtered alignments)
INPUT_SUBFOLDER="filtered_50%Taxa_50Col_50%ColumnGapTolerance_50%RowGapTolerance"

# Set working directories
INPUT_DIR="$PROJECT/input/$INPUT_SUBFOLDER"
# Output subfolder named after input filter
OUTPUT_DIR="$PROJECT/output/$INPUT_SUBFOLDER"
# Script and AMAS location
SCRIPT_DIR="$PROJECT/scripts"
AMAS_SCRIPT="$SCRIPT_DIR/AMAS.py"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Move to the directory containing AMAS.py
cd "$(dirname "$0")" || exit

# Run AMAS concatenation
python "$AMAS_SCRIPT" concat \
  -i "$INPUT_DIR"/*.fasta \
  -f fasta \
  -d dna \
  --concat-part "$OUTPUT_DIR/partitions.txt" \
  --concat-out "$OUTPUT_DIR/concatenated.fasta" \



# Convert to nexus
python "$AMAS_SCRIPT" convert \
  -d dna \
  -f fasta \
  -i "$OUTPUT_DIR/concatenated.fasta" \
  -u nexus


# Move the concatenated nexus file to the output directory
mv "$SCRIPT_DIR/concatenated.fasta-out.nex" "$OUTPUT_DIR/"