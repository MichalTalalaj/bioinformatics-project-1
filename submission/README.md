# Bioinformatics - Project 1

Use Python 3.12. From this folder, install dependencies:

```sh
python -m pip install -r requirements.txt
```

## Assignment 1

```sh
cd assignment_1
python orf_finder.py genome.fasta --min-len 100 --out orfs.csv
python kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv
cd ..
```

## Assignment 2

Requires internet access. The generated example report is included.

```sh
python assignment_2/genereport.py APOE --organism "Homo sapiens" --out assignment_2/reports/
```

## Assignment 3

```sh
cd assignment_3
python -m notebook rnaseq.ipynb
```

Select the installed Python environment as the kernel and run all cells from the beginning. All input files are included; keep the notebook's working directory as assignment_3/. Do not rerun the log-transformation cell without reloading the input.

## Analysis limitations

- Driver-mutation metadata is unavailable: mutation-coloured plots and EGFR/KRAS tests are omitted; smoking-cluster tests are used.
- The 77 tumour-normal candidate pairs have matching available clinical metadata; patient identities were not independently verified.