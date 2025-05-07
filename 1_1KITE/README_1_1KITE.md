# Step 1: 1KITE Phylogenomic Pipeline

This directory contains all scripts, software, and outputs for the **1KITE pipeline**, adapted for processing hybrid-capture Illumina reads from *Romalea* species.

The pipeline performs:
- Adapter trimming and quality filtering
- De novo gene assembly
- Orthology inference using Orthograph
- Alignment and cleaning of orthologous loci

---

## Folder Structure

```
1_1KITE/
├── Output/              # Final FASTA alignments of individual orthologous genes (for downstream phylogenetic analysis)
│
├── Software/            # Core software used in the pipeline
│   ├── trimmomatic/         # Adapter removal and quality filtering
│   ├── SOAPdenovo2/         # De novo transcriptome assembly of cleaned reads
│   └── Orthograph/          # Orthology inference via HMM-based reciprocal search
│
├── scripts_1kite/       # All custom scripts and intermediate results
│                        # Includes wrapper scripts for gene selection, alignment trimming, and formatting
│
├── README_AHE__data_analysis_notes_upstream_TAMU_Edition_updated.sh
│                        # Master workflow guide with detailed step-by-step instructions for running the pipeline
```

---

## Software Tools

- **Trimmomatic**  
  Removes adapter sequences and performs quality filtering on raw Illumina reads.

- **SOAPdenovo2**  
  Assembles the filtered reads de novo into contigs or transcriptomes.

- **Orthograph** ([GitHub](https://github.com/mptrsen/Orthograph))  
  Performs graph-based orthology inference using profile HMMs against a reference probe set (in this case, Orthoptera-specific markers).

---

## Notes

- The pipeline is adapted from the *1KITE* phylogenomic protocol, widely used for transcriptome-based phylogenetic inference in insects.
- Ortholog detection is tailored to an **Orthoptera-specific probe set** for higher recovery and accuracy.
- Refer to the shell script `README_AHE__data_analysis_notes_upstream_TAMU_Edition_updated.sh` for the full step-by-step execution workflow, expected input/output file formats, and command structure.

---

## Important

- Make sure to configure paths and environment variables inside the scripts depending on your local or cluster environment.
- Intermediate outputs are saved throughout to allow for inspection, troubleshooting, or rerunning individual steps.

---
