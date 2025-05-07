# Romalea Phylogenomic Pipeline

This repository contains pipeline for performing phylogenomic inference on species of the genus *Romalea*. The pipeline uses Anchored Hybrid Enrichment data to construct a species tree, and it includes steps for read processing(Trimmomatic), orthology inference (Orthograph), alignment (MAFFT), filtering (filter_alignment.py), concatenation (AMAS.py), and maximum likelihood inference using **IQ-TREE v2**.

---

## Project Structure

```
Romalea_phylogenomics/
├── 1_1KITE/       # Raw read QC, de novo assembly, orthology inference, gene alignment
├── 2_filter/      # Filters alignments by length, taxon completeness, and missing data
├── 3_AMAS/        # Concatenates alignments and generates partition file using AMAS
├── 4_IQ_TREE/     # Runs IQ-TREE on concatenated alignments for ML inference
```

Each subfolder contains its own `README.md` with step-specific instructions.

---

## Pipeline Overview

| Step | Description |
|------|-------------|
| `1_1KITE/` | Performs quality control (Trimmomatic), de novo assembly (SOAPdenovo2), ortholog prediction (Orthograph), and gene alignment (MAFFT). Outputs individual gene alignments in FASTA format. |
| `2_filter/` | Filters alignments to remove short genes, poorly sampled taxa, and gap-rich columns using custom Python script (filter_alignment.py). |
| `3_AMAS/` | Concatenates filtered alignments using [AMAS](https://github.com/marekborowiec/AMAS) and creates a partition file for downstream analyses. |
| `4_IQ_TREE/` | Performs model selection and tree inference with ultrafast bootstrap using IQ-TREE v2. Uses NEXUS alignment and partition file as input. |

---

## Software Requirements

- **Trimmomatic** (read trimming)
- **SOAPdenovo2** (de novo assembly)
- **Orthograph** ([link](https://github.com/mptrsen/Orthograph)) (orthology inference)
- **AMAS** ([link](https://github.com/marekborowiec/AMAS)) (concatenation and partitioning)
- **IQ-TREE v2** (phylogenetic inference)
- **Python 3 + Biopython** (for filtering scripts)
- **Conda** (recommended for environment management)

---

## General Usage

The pipeline is modular, and each step must be run sequentially:

1. **1_1KITE/**: Assemble and align orthologous genes from raw Illumina reads.
2. **2_filter/**: Clean alignments by filtering by row/column completeness, length, and gap content.
3. **3_AMAS/**: Concatenate alignments and generate partition files using AMAS.
4. **4_IQ_TREE/**: Run ML phylogenetic inference using IQ-TREE v2 on the supermatrix.

> Note: The pipeline is **not fully automated**. Users must edit file paths and parameters within each step's script(s).

---

## Output Summary

| Step | Outputs |
|------|---------|
| `1_1KITE/` | FASTA alignments of individual orthologous genes |
| `2_filter/` | Filtered FASTA gene alignments |
| `3_AMAS/` | Concatenated alignment in FASTA and NEXUS format and partition file |
| `4_IQ_TREE/` | ML tree files, bootstraps, and model reports |

---

## Reproducibility

In 2_filter/ and onwards, the scripts logs filtering parameters and outputs using a consistent folder naming scheme (e.g., `filtered_50%Taxa_50Col_50%Gap`). This allows re-running or comparing filtered datasets and tree topologies.

A future goal is to wrap this pipeline in a automatized workflow to improve reproducibility and scalability on high-performance computing clusters.

---

## Contact

For questions or collaboration, please contact hsong@tamu.edu, jorgemedinad@tamu.edu or jackson.linde@tamu.edu
