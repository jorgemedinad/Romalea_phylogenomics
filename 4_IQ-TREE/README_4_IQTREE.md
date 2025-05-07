# IQ-TREE Phylogenomic Pipeline

This step contains scripts and a standardized folder structure to automate phylogenomic inference using **IQ-TREE v2**. The pipeline performs model selection and tree inference with ultrafast bootstrap support on concatenated sequence alignments and partition files.

---

## Directory Structure

```
4_IQ-TREE/
├── input/
│   └── filtered_<params>/                # Contains filtered alignments and partition files
│       ├── concatenated.fasta-out.nex    # Concatenated alignment (NEXUS format)
│       └── partitions.txt                # Partition file with ranges per gene
│
├── output/
│   └── filtered_<params>/                # Output directory corresponding to each input set
│       ├── model_selection.*             # Model selection and merged partition outputs
│       └── bootstrap_tree.*              # Tree inference with bootstrap support
│
├── scripts/
│   └── run_iqtree.sh                     # Bash script to run the full pipeline
```

---

## Requirements

- IQ-TREE v2
- Miniconda (recommended)
- Bash

Create and activate a conda environment with IQ-TREE:
```bash
conda create -n iqtree_env -c bioconda iqtree
conda activate iqtree_env
```

---

## Running the Pipeline

### 1. Activate your Conda environment
```bash
conda activate iqtree_env
```

### 2. Run the script
```bash
bash scripts/run_iqtree.sh <filtered_input_folder_name>
```

#### Example:
```bash
bash scripts/run_iqtree.sh filtered_50%Taxa_50Col_50%ColumnGapTolerance_50%RowGapTolerance
```

This will:
- Prepend `"DNA, "` to each line in `partitions.txt`
- Run **ModelFinder** with partition merging (`MF+MERGE`) using "iqtree2 -s concatenated_alignment.nex -p partitions.txt -m MF+MERGE -T AUTO"
- Run **IQ-TREE** with ultrafast bootstrap (1000 replicates) using "iqtree2 -s concatenated_alignment.nex -p model_selection.best_scheme.nex -B 1000 -T AUTO"
- Save outputs to:  
  `output/<filtered_input_folder_name>/`

---

## Notes

- Alignment must be in NEXUS format.
- Partition file should use simple numeric ranges per gene.
- All output files are named using the input folder name for reproducibility.
