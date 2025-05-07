#!/bin/bash

# === Bash Wrapper for Filtering nexus concatenated fileScript ===
# Define filtering parameters
TOTAL_TAXA=33 # write the  number of taxa of your specific analysis
COLUMN_GAP=0.50 # maximum proportion of gaps allowed per column (larger percentages are less stringent)
ROW_GAP=0.50 # maximum proportion of gaps allowed per row (i.e., taxa) (larger percentages are less stringent)
MIN_COLUMNS=50 # minimum length of alignment (minimum number of columns in alignment matrix)
TAXA_PERCENT=0.5 # minimum proportion of taxa to retain in alignments

# Locate the script relative to this wrapper's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/filter_alignments.py"

# Execute the Python script with arguments
python3 "$PYTHON_SCRIPT" \
  --total_taxa $TOTAL_TAXA \
  --column_gap $COLUMN_GAP \
  --row_gap $ROW_GAP \
  --min_columns $MIN_COLUMNS \
  --taxa_percent $TAXA_PERCENT
