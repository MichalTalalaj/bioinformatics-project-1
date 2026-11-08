# Bioinformatics — Project 1

Tested with Python 3.12. Install dependencies from the project root using your chosen Python environment:

```sh
python -m pip install -r requirements.txt
```

## Assignment 1 — ORFs and k-mers

Run the assignment's example commands from its directory:

```sh
cd assignment_1
python orf_finder.py genome.fasta --min-len 100 --out orfs.csv
python kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv
```

The included FASTA is lambda bacteriophage `NC_001416.1` (48,502 nt). Both commands write or replace the specified output CSV.

- `orfs.csv`: frame, 1-based inclusive positions on the original forward reference, nucleotide length including STOP, and protein without the terminal `*`. The minimum length is in nucleotides.
- `kmers_k3.csv`: overlapping forward-sequence 3-mers, counts and relative frequencies, sorted by decreasing count with alphabetical ties.

The scripts accept one nonempty FASTA record, convert DNA to uppercase and reject symbols outside `ACGT`. Run either script with `--help` for its arguments. Keep both scripts together: `kmer_freq.py` reuses the FASTA loader in `orf_finder.py`.

## Assignment 2 — gene report

From the project root, enter the assignment directory and run:

```sh
cd assignment_2
python genereport.py APOE --organism "Homo sapiens" --out reports/
```

Internet access is required. The program retrieves NCBI Gene, linked RefSeq mRNA, UniProtKB and Ensembl annotations, then writes `reports/APOE_Homo_sapiens.md`. It creates the output directory if needed and replaces the same report on successful reruns. A failed lookup leaves an existing report unchanged.

The report includes all linked RefSeq mRNA variants and prefers all exact reviewed UniProt matches, falling back to unreviewed matches when needed. Missing annotations are labelled `Unavailable`. CDS locations are transcript-relative; Ensembl genomic coordinates include their assembly and strand. Both use 1-based inclusive positions. Run `python genereport.py --help` for the arguments.

## Assignment 3

Not implemented yet; its directory is reserved for later work.

## Reference

The assignment description is in [docs/Bioinformatics_Project_1.pdf](docs/Bioinformatics_Project_1.pdf). Planning and study notes are also in `docs/`; `AGENTS.md` and `requirements.txt` remain at the project root.
