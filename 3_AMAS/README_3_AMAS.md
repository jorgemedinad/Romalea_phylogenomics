# Running AMAS Concatenation Script

This step uses [AMAS](https://github.com/marekborowiec/AMAS) to concatenate multiple gene alignments in FASTA format into a single concatenated alignment with an associated partition file.

## Directory Structure

```
3_AMAS/
  ├── input/                      # Contains input filtered gene alignments (*.fasta)
  ├── output/                     # Output folder for concatenated files
  │ └── filtered_<params>/        # Concatenated files stored in folder named by parameter set
  ├── scripts/                    # Contains AMAS.py and this run script
    └── run_amas.sh               # Bash script to execute AMAS
    └── AMAS.py                   # AMAS script 
```

## How to Run

1. Make sure AMAS and Python are properly installed in your environment.
2. Place all individual filtered gene FASTA alignments inside the `input/` directory.
3. Open a terminal and run the script:

```bash
bash scripts/run_amas.sh
```

## What the Script Does

- Runs `AMAS.py` using the `concat` function.
- Searches for all `*.fasta` files in the `input/` directory.
- Concatenates all alignments assuming **DNA data** in **FASTA** and **NEXUS** format.
- Produces the following output files in the `output/` directory:
  - `concatenated.nex` – concatenated alignment in **FASTA** format.
  - `partitions.txt` – partition file in **RAxML** format.
  - `concatenated.fasta-out.nex` – partition file in **NEXUS** format.

## Dependencies

- Python 3
- [AMAS](https://github.com/marekborowiec/AMAS) (should be located in `scripts/AMAS.py`)
- Make sure Python can run `AMAS.py` with:

```bash
python scripts/AMAS.py --help
```
