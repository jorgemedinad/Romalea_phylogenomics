# Running Alignment Filtering Script

This script filters gene alignments in FASTA format obtained from the 1KITE pipeline to improve data quality by removing poorly sampled positions and taxa with excessive missing data. It ensures alignments meet defined thresholds before concatenation with [AMAS](https://github.com/marekborowiec/AMAS).

---

## Directory Structure

```
2_Filter/
├── input/                               # Raw alignments from 1KITE pipeline (one per gene, FASTA format)
├── output/
│   └── filtered_<params>/               # Filtered alignments named by parameter set
├── scripts/
│   ├── run_filtering.sh                 # Bash wrapper for batch filtering
│   └── filter_alignment.py              # Python script performing the actual filtering
```

---

## How to Run

1. Place all unfiltered alignments in the `input/` directory.
2. Open and modify `scripts/run_filtering.sh` to set your desired thresholds.
3. Then run the script:

```bash
bash scripts/run_filtering.sh
```

Filtered alignments will be saved to a subfolder inside `output/`, named according to the parameter set used.

---

## Example Output

If using 50% minimum taxon coverage per alignment, 50% maximum column gap, and 50% maximum row gap per taxon, output will go to:

```
output/filtered_50%Taxa_50Col_50%ColumnGapTolerance_50%RowGapTolerance/
```

Each filtered file retains the original gene name.

---

## What the Script Does

For each alignment:
- **Removes taxa (rows)** that exceed a user-defined **row gap tolerance** (i.e., taxa with too many gaps/missing data).
- **Removes positions (columns)** where the number of taxa with data is less than a defined **minimum taxon threshold** or the percentage of gaps exceeds the **column gap tolerance**.
- **Skips alignments** that are shorter than a minimum number of positions or contain too few taxa.
- Writes a filtered version of each alignment in FASTA format.

Only alignments that pass these filters are retained.

---


## Parameters Explained

Set these in `run_filtering.sh`:

| Parameter        | Description |
|------------------|-------------|
| `TOTAL_TAXA`     | Total number of taxa expected in each alignment. Used to compute proportions for filtering. |
| `COLUMN_GAP`     | Maximum proportion of gaps allowed per **column** (alignment position). For example, `0.50` allows columns where up to 50% of taxa have gaps. |
| `ROW_GAP`        | Maximum proportion of gaps allowed per **taxon**. Taxa with more missing data than this threshold are removed. |
| `MIN_COLUMNS`    | Minimum number of alignment columns (i.e., minimum alignment length) required to retain a gene after filtering. |
| `TAXA_PERCENT`   | Minimum proportion of the total taxa (`TOTAL_TAXA`) that must remain in the alignment after row and column filtering. If too many taxa are removed, the alignment is skipped. |

---

## Example Configuration

In `scripts/run_filtering.sh`:

```bash
TOTAL_TAXA=33
COLUMN_GAP=0.50
ROW_GAP=0.50
MIN_COLUMNS=50
TAXA_PERCENT=0.5
```

This keeps alignments where:
- Each site (column) has data for at least 50% of taxa.
- Each taxon (row) has at least 50% non-gap characters.
- The resulting alignment has at least 50 taxa and 50 columns.

---

## Virtual Environment

This script is intended to be run within a Python virtual environment to ensure clean dependency management. For example:

```bash
conda create -n amas_filter_env python=3.10
conda activate amas_filter_env
pip install biopython
```

Activate the environment before running the filtering script:

```bash
conda activate amas_filter_env
bash scripts/run_filtering.sh
```

This ensures compatibility with Biopython and isolates dependencies from your system Python installation.

---

## Next Step: Concatenation

Once filtering is complete, use the filtered output as input for `AMAS.py` in the concatenation pipeline.
