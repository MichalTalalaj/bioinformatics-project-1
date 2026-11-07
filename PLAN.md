# Bioinformatics 2026 — Project 1: staged implementation plan

**Authoritative requirements:** `Bioinformatics_Project_1.pdf` (3 pages).  
**Workflow:** approval required **before coding each stage**, followed by a second approval **before proceeding**. Codex must stop between stages.  
**Priority:** meet every PDF requirement with the simplest correct, understandable code. Do not add optional functionality.

## Initial decisions (before coding)

- [x] **S0 — Inspect inputs and agree on conventions.** Confirm that the PDF is available in the working directory. Confirm Python environment and only necessary dependencies. Choose a reference genome (default proposal: lambda bacteriophage `NC_001416.1`). Decide and record minimal conventions where the assignment is silent: ORF length units and translation terminator, handling ambiguous symbols, output format for gene report, output structure for RNA-seq, basic statistical defaults. Do not implement anything at this stage.

In every stage below: **propose → obtain approval → implement only that stage → verify → explain → await acceptance**. Mark a checkbox only after acceptance.

### S0 — approved conventions and inspection findings

The student accepted S0 and separately approved environment setup and A1.1. These are project conventions where the PDF leaves a choice, not additional assignment requirements. Later stages still require their own approval.

- **Inputs:** `Bioinformatics_Project_1.pdf` is present, and all three pages have been read alongside `AGENTS.md` and this plan.
- **Environment:** Python 3.12.14 runs at `D:\Anaconda\envs\ds-semester\python.exe`. Module discovery found `requests`, NumPy, pandas, SciPy, Matplotlib and scikit-learn; it did not find Biopython (`Bio`) or `umap`. This confirms discovery only, not successful imports or full functionality. No dependencies were installed. Biopython is required for Assignment 1; UMAP support is needed for Assignment 3. Installation requires separate approval.
- **Reference genome:** lambda bacteriophage `NC_001416.1`. Download and inspect the FASTA in A1.1.
- **ORF length:** count nucleotides, including both the initial `ATG` and terminal STOP codon. Apply `--min-len` in these same units.
- **ORF starts and stops:** consider every `ATG` as a separate candidate, ending at its first in-frame `TAA`, `TAG` or `TGA`. Nested ORFs are allowed. Do not report candidates without a terminal STOP.
- **Coordinates and frames:** use 1-based, inclusive positions on the original forward reference, with `start <= end` on both strands. Use frames `+1`, `+2`, `+3`, `-1`, `-2`, `-3`; negative frames refer to offsets in the reverse complement. The frame identifies orientation.
- **Translation:** use the standard genetic code and omit the terminal `*` from the protein sequence because STOP is not an amino acid. The reported nucleotide length still includes the STOP codon.
- **Ambiguous symbols in Assignment 1:** convert input sequence letters to uppercase, then reject any symbol outside `ACGT` with a clear error. Do not guess bases or silently discard sequence windows.
- **k-mers:** count overlapping windows in the supplied forward sequence, without additionally counting its reverse complement. Report observed k-mers. For a sequence of length `N` and valid `1 <= k <= N`, relative frequency is `count / (N - k + 1)`. Sort by decreasing count, breaking ties alphabetically.
- **Gene report:** use `APOE` and `Homo sapiens` as the example. Produce one Markdown report in the directory selected by `--out`, marking absent fields as unavailable. Agree on record-selection rules in the relevant stages after inspecting real API responses.
- **RNA-seq output structure:** use one `rnaseq.py` script, extended only in approved stages, and a `results/` directory for required tables and plots.
- **Statistical defaults:** use `random_state=42` where supported for stochastic methods. For differential expression, use Benjamini-Hochberg FDR correction with a significance threshold of 0.05. Specific tests, expression filtering, method parameters and cross-dataset compatibility remain subject to inspection and approval in A3.0.
- **S0 scope:** record decisions only in `PLAN.md`; do not implement scripts, download datasets or install packages. Add study-guide explanations alongside approved implementations.

### A1.1 — implementation status (accepted)

- Created `.venv/` using Python 3.12.14. The student installed Biopython manually; version 1.88 was confirmed in that environment and pinned in `requirements.txt` as the only direct dependency. `pip check` reported no broken requirements.
- Added `load_genome()` in `orf_finder.py` to read exactly one nonempty FASTA record, uppercase its sequence and reject symbols outside `ACGT`.
- After command-environment network failures, the student supplied `genome.fasta`. The loader read one record with ID `NC_001416.1`, description `Enterobacteria phage lambda, complete genome`, length 48,502 nucleotides and only `ACGT` symbols. This verifies the local input's header and content structure; it is not an independent checksum comparison with NCBI.
- Functional checks passed for the real genome, mixed-case multiline DNA, and rejection of missing files, no records, empty sequences, multiple records, ambiguous bases, invalid symbols and gaps. Small test sequences were supplied in memory; the genome was left unchanged.
- Python syntax and Git whitespace checks passed; `.venv/Scripts/python.exe` is ignored. The student accepted A1.1 and approved the proposed A1.2 implementation.

### A1.2 — implementation status (accepted)

- Added `find_orfs(sequence)` in `orf_finder.py`. It checks every `ATG` in all six frames and reports the first in-frame STOP, signed frame, 1-based inclusive original-reference coordinates, nucleotide length including STOP and protein without the terminal `*`. Biopython performs reverse complement and translation only; the ORF search is implemented explicitly.
- Hand-checked cases passed for all six frames, reverse coordinates, nested starts, all three STOP codons, first-STOP behavior, out-of-frame stops, sequence boundaries, separate ORFs, trailing incomplete codons and missing complete ORFs.
- The lambda genome produced 1,800 ORFs without minimum-length filtering: `+1`: 340, `+2`: 332, `+3`: 326, `-1`: 283, `-2`: 248, `-3`: 271. Every result passed coordinate-to-sequence reconstruction, frame, length, ATG, first-STOP and protein checks; no duplicate frame/coordinate entries were found. These are sequence-defined ORFs, not confirmed genes.
- Added the method, coordinate formulas, verification results and instructor questions to `STUDY_GUIDE.md`. The student accepted A1.2. CLI arguments, minimum-length filtering and CSV output are covered in A1.3 below.

### A1.3 — implementation status (accepted)

- Added a `main()` entry point to `orf_finder.py` using standard-library `argparse` and `csv`. It accepts a positional FASTA path, `--min-len` (positive integer, default 100 nucleotides including STOP) and `--out` (default `orfs.csv`). It retains ORFs with `length >= min_len` and writes `frame,start,end,length,protein` columns. Existing output files are preserved; likely input/output errors produce short messages and exit code 2.
- Ran `.venv\Scripts\python.exe -B orf_finder.py genome.fasta --min-len 100 --out orfs.csv`. The output contains 708 ORFs, with lengths from 102 to 3,399 nucleotides. Counts by frame: `+1`: 146, `+2`: 129, `+3`: 138, `-1`: 110, `-2`: 92, `-3`: 93.
- CLI checks passed for inclusive threshold boundaries, 100-nt/default-output behavior, a header-only result, reverse-strand CSV fields, lowercase input, invalid thresholds, missing/invalid FASTA, a missing output directory, existing-file protection and `--help`.
- Every generated CSV row matched the filtered search results and passed coordinate, frame, length, ATG, first-STOP and translation checks. The student accepted A1.3. The k-mer implementation is covered in A1.4 below.

### A1.4 — implementation status (accepted)

- Added `kmer_freq.py`, reusing the Biopython FASTA loader from `orf_finder.py`. `count_kmers(sequence, k)` uses standard-library `Counter` to count overlapping windows of the forward sequence and returns `(kmer, count, frequency)` rows. It requires `1 <= k <= N`, uses `count / (N - k + 1)` for frequencies and sorts by decreasing count with alphabetical ties.
- The CLI accepts a positional FASTA path, required integer `-k` and required `--out`. It writes `kmer,count,frequency` columns, preserves existing output files and gives short messages for likely errors. No additional dependencies were introduced.
- Ran `.venv\Scripts\python.exe -B kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv`. The output contains 64 distinct 3-mers with a count sum of 48,500 (`48,502 - 3 + 1`); relative frequencies sum to 1 within floating-point tolerance. All counts were independently checked against overlapping occurrences in the genome, and every row's frequency and sort position were verified.
- Hand-counted examples passed for overlapping windows, unequal counts, alphabetical ties, `k=1`, `k=N`, forward-only counting and lowercase multiline input. CLI argument, invalid/missing FASTA, output-path, existing-file and help checks passed. At the student's request, counting and row construction were expanded into explicit loops with a named sorting rule; behavior-preservation checks passed. The student accepted A1.4. The Assignment 1 audit is covered in A1.5 below.

### A1.5 — acceptance check (verified; awaiting acceptance)

- Re-read Assignment 1 directly from page 1 of `Bioinformatics_Project_1.pdf` and compared every requirement with the scripts, conventions and outputs. No implementation gaps or necessary corrections were identified.
- Confirmed all four required deliverables exist and are nonempty: `orf_finder.py`, `kmer_freq.py`, `orfs.csv` and `kmers_k3.csv`. Both scripts use the same Biopython loader and the local lambda input identified as `NC_001416.1`. The ORF search is implemented explicitly; Biopython is used for parsing, reverse complement and translation.
- Ran both assignment-style CLIs with the project interpreter into a temporary directory. Each regenerated CSV was byte-for-byte identical to the corresponding saved CSV. SHA-256 checks confirmed that the scripts, input genome and saved tables were unchanged by the audit; temporary outputs were cleaned up.
- Rechecked exact examples for all six frames and all 708 saved ORF rows: required fields, minimum length, 1-based inclusive coordinates, frame mapping, ATG, first in-frame STOP and translation. Independently recounted every one of the 64 saved 3-mers, confirming 48,500 windows, frequencies and sorting. Both CLI help commands and `pip check` passed.
- Recorded the requirement mapping and results in `STUDY_GUIDE.md`. No source code or saved CSV changes were needed. A1.5 awaits acceptance; its checkbox remains unchecked. A2.1 has not started.

## Assignment 1 — ORFs and k-mer frequencies (5 points)

**Explicit deliverables:** `orf_finder.py`, `kmer_freq.py`, `orfs.csv`, `kmers_k3.csv` (plus an input FASTA genome for local runs).

- [x] **A1.1 — Input genome and FASTA reading.** Obtain the chosen NCBI genome and inspect its FASTA record; implement the minimum Biopython FASTA loading needed by the ORF script. Verify input is read correctly. Keep files minimal.
- [x] **A1.2 — ORF search.** Implement your own ATG → first in-frame stop search in all 3 forward and 3 reverse-complement frames (no Biopython ORF-finding function). Calculate 1-based original-sequence positions, lengths, frame, and translated proteins. Verify on tiny hand-checked forward and reverse examples, including boundaries.
- [x] **A1.3 — ORF CLI and CSV.** Add required `--min-len` and `--out` interface; write `orfs.csv` for the genome. Validate columns and a sample of coordinates and translations. Example: `python orf_finder.py genome.fasta --min-len 100 --out orfs.csv`.
- [x] **A1.4 — k-mer script and CSV.** Read FASTA with Biopython, accept `-k` and `--out`, count overlapping sequence k-mers, compute relative frequencies, order by descending count; generate `kmers_k3.csv`. Verify on an easily counted toy DNA string. Example: `python kmer_freq.py genome.fasta -k 3 --out kmers_k3.csv`.
- [ ] **A1.5 — Assignment 1 acceptance check.** Confirm both scripts run and both requested CSVs exist, are nonempty where expected, and conform to the PDF. Review only necessary corrections.

## Assignment 2 — automated gene report via API (5 points)

**Explicit deliverable:** `genereport.py`, plus the generated report for the selected gene/organism. The PDF does not mandate a specific report file type; agree on **one simple format** (suggested: Markdown or text), not multiple formats.

- [ ] **A2.1 — Minimal interface and NCBI Gene.** Agree on gene/organism example (PDF: APOE, Homo sapiens). Implement CLI and `requests` access to NCBI `esearch.fcgi` (`db=gene`) and `esummary.fcgi`. Extract Gene ID, full name, chromosomal location, aliases, summary, MIM references. Inspect real response to verify field mapping.
- [ ] **A2.2 — NCBI Nucleotide.** Use `elink.fcgi` (gene → nuccore) and `efetch.fcgi` (GenBank / FASTA as needed), using `requests` only. Extract RefSeq mRNA accession/version, length, last modification date, CDS coordinates, `protein_id`, `transl_table`. Confirm correct record choice and parse only the necessary fields, without Biopython.
- [ ] **A2.3 — UniProtKB.** Query `rest.uniprot.org/uniprotkb/search` (`format=json`, appropriate `fields`) using `requests`. Include accession, reviewed status, protein name, function, length, ligands/cofactors, diseases, PDB structures, keywords. Verify gene/organism match and handle legitimately missing information.
- [ ] **A2.4 — Ensembl.** Query `rest.ensembl.org/lookup/symbol/{species}/{symbol}?expand=1` via `requests`. Extract gene ID, coordinates, transcript count, canonical transcript, biotype. Verify organism name formatting and returned data.
- [ ] **A2.5 — Assemble and verify report.** Produce one readable report containing every requested field, use `--out` directory, run `python genereport.py APOE --organism "Homo sapiens" --out reports/`. Check each database section against API results. Do not invent missing values.

## Assignment 3 — lung adenocarcinoma RNA-seq (10 points; 2 per subtask)

**Required sources:** GSE40419 (RPKM supplementary matrix and sample metadata); GSE10072 for comparisons and external classification testing.  
**Outputs:** a small number of readable analysis files/plots/results sufficient to demonstrate each requested result. The PDF does **not** prescribe script names or figure filenames for this assignment; decide these before implementation and avoid creating an unnecessary file per operation.

- [ ] **A3.0 — Inspect datasets and analysis choices (no analysis code yet).** Download/inspect GSE40419 RPKM and metadata and GSE10072 from their authentic source files; match sample and patient IDs; identify clinical variable coverage (`tissue`, `patient`, `sex`, `smoking`, `driver mutation`, including EGFR/KRAS); identify tumour–normal pairs. Check how expression features align across datasets. Propose minimal statistical/technical conventions, and explicitly flag missing metadata or incompatibilities to the student before continuing.
- [ ] **A3.1 — Preprocessing and sample quality (2 pts).** Load matrix and metadata; apply `log2(RPKM + 1)`; filter low-expression genes using an explicitly stated, approved simple rule; produce a sample–sample correlation heatmap. Verify samples/genes, mapping, and transformations.
- [ ] **A3.2 — Dimensionality reduction (2 pts).** Generate PCA, t-SNE, and UMAP sample plots coloured by tissue **and** available clinical variables; compare them with silhouette scores. Use reproducible settings where applicable. Explain what each plot/score actually supports.
- [ ] **A3.3 — Clustering (2 pts).** Apply k-means and hierarchical clustering, compute ARI against tissue labels; repeat clustering on tumour samples only; test association of tumour clusters with smoking and EGFR/KRAS status using a suitable, stated statistical test. Report counts and results clearly, including missing metadata if any.
- [ ] **A3.4 — Paired differential expression (2 pts).** Construct verified tumour–normal patient pairs; run an appropriate paired gene-wise test; apply FDR correction; create volcano plot; compare results with GSE10072 from Lab 3. Confirm comparability before interpreting any overlap/differences.
- [ ] **A3.5 — Classification and external test (2 pts).** Build a **simple** tumour-vs-normal classifier. Use patient-grouped cross-validation with strictly training-only learned preprocessing and feature selection. Report cross-validated ROC AUC and most important genes (with a clear definition of importance). Test the final model on GSE10072 only after confirming a compatible feature/scale representation; report external results and limitations. Never fit or tune using GSE10072 labels.
- [ ] **A3.6 — Assignment 3 acceptance check.** Ensure that all five parts have their requested tables/plots/statistics and that no claim relies on missing data, mismatched samples, or leakage. Fix only gaps against the PDF.

## Final review and study preparation

- [ ] **F1 — Requirement-by-requirement audit.** Compare the implementation and outputs with each line of the PDF. List any gaps, with no optional extensions or extra features.
- [ ] **F2 — Complete `STUDY_GUIDE.md`.** For each *implemented* part: inputs/outputs, step-by-step explanation of the code and method, interpretation of results, verification, and short answers to likely instructor questions. The guide is **study material**, not a mandatory submission item.
- [ ] **F3 — Submission inventory.** Identify exactly the files required by the PDF and any additional files genuinely needed to present Assignment 2–3 results. Do not assume that the PDF requires `AGENTS.md`, `PLAN.md`, `STUDY_GUIDE.md`, a README, or unrelated infrastructure to be submitted.

## Known ambiguities — do not silently decide

The PDF does not specify all conventions. Before relevant implementation, agree on a minimal defensible choice for:
- ORF minimum-length units, multiple possible start codons within a region, translation stop marker, and ambiguous-base handling.
- Relative-frequency denominator/ambiguous k-mers (use a clear, consistent rule).
- Report output type and response selection among multiple gene, RNA, protein, and transcript matches.
- Exact GSE40419 metadata mapping, expression filtering threshold, t-SNE/UMAP settings, cluster counts, association test choice, classifier/importance method.
- GSE10072 comparison method, cross-platform gene identifiers, compatible input scale, and what to do if a requested field is not available.

**Rule:** Make assumptions explicit, prefer the simplest correct option, and obtain the student's approval at the corresponding stage.
